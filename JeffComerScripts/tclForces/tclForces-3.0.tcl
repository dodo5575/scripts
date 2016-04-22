## this tclforces script defines a set of procedures that will be called in external scripts each calcforces cycle

## This code might be faster by creating namespaces for each proc in each forceFile sourced by it.  Then the namespace can store information that doesn't need to be recalculated at every timestep (such as total mass).

## change this for your filesystem
set tclForcesDir $env(HOME)/scripts/tclForces
set forceDefFiles "bonds.1.tcl restraints.1.tcl forces-2.1.tcl smdForces.1.tcl DNAforces.1.tcl"

source $tclForcesDir/vector.tcl

namespace eval ::tclForces {
    ## Namespace that contains procs/vars that interface between (relatively) abstract, high-level force applying procs and the low-level interface provided by namd for actually applying the forces
    
    ## import/export
    namespace import vector::*
    namespace export getCOM getCOMZ getMass sumWeight addEnergy applyForce printForce applyTorque applyTorques
    
    ## define arrays and variables
    array set forces ""
    array set energy ""
    array set coords ""
    array set masses ""

    variable energyNames ""    
    variable forceFileContents ""
    variable atoms ""
    variable forcecount 2
    variable printCount 0
    variable timeStep 0
    if {[info exists ::firstTimeStep]} { set timeStep $::firstTimeStep }

    #### define public procs to get atom positions
    proc getCOM { ids } {
	foreach var "coords masses" { variable $var }
	if {[llength $ids] == 1} { return "{$coords($ids)} $masses($ids)" }
	set c ""
	set m 0
	foreach id $ids {
	    lappend c [vecscale $coords($id) $masses($id)];# build a list of {m1*(x1,y1,z1) ,  ...} and add elements later to get COM
	    set m [expr {$m+$masses($id)}]
	}
	set c [vecscale [eval vecadd $c] [expr {1./$m}]] ;# COM!
	return "{$c} $m"
    }
    proc getCOW { ids {weightName masses} } {
	variable coords
	variable $weightName
	set c ""
	set w 0
	foreach id $ids {
	    set weight [subst $${weightName}($id)]
	    lappend c [vecscale $coords($id) $weight];# build a list of {w1*(x1,y1,z1) ,  ...} and add elements later to get COW
	    set w [expr {$w+$weight}]
	}
	set c [vecscale [eval vecadd $c] [expr {1./$w}]] ;# COW!
	return "{$c} $m"
    }
    proc getCOMZ { ids } {
	foreach var "coords masses" { variable $var }
	if {[llength $ids] == 1} { return "[lindex $coords($ids) 2] $masses($ids)" }
	set c ""
	set m 0
	foreach id $ids {
	    lappend z "[lindex $coords($id) 2]*$masses($id)"
	    set m [expr {$m+$masses($id)}]
	}
	set z [expr [join $z "+"]]
	set z [expr {$z/$m}] ;# COM!
	return "$z $m"
    }
    proc getMass { ids } {
	foreach var "masses" { variable $var }
	set m 0
	foreach id $ids { set m {$m+$masses($id)} }
	return $m
    }
    proc sumWeight { ids weightName } {
	variable $weightName
	set sum 0
	foreach id $ids {
	    append sum [subst "+$${weightName}($id)"]
	}
	set sum [expr $sum]
    }
    ## I'm not sure that the next two procs are any good...
    proc getWeights { idName ids code } {
	upvar $idName id
	set weights ""
	foreach id $ids {
	    lappend weights [uplevel $code]
	}
	return $weights
    }
    proc sumWeight2 { idName ids code } {
	upvar $idName id
	set sum 0
	foreach id $ids {
	    append sum "+[uplevel expr $code]"
	}
	return [expr $sum]
    }
    ## procs to apply forces/maintain a list of energy and print force information
    proc applyForce {ids force {weightName masses} {weightSum ""} {print 0}} {
	variable forces
	variable $weightName
	if {$weightSum == ""} { set weightSum [sumWeight $ids $weightName] }
	
	set force [vecscale [expr {1./$weightSum}] $force]
	foreach id $ids {
	    set forces($id) [vecadd $forces($id) [vecscale [subst $${weightName}($id)] $force]]
	}
    
    }
    proc applyEachTorque {ids torque {axis "0 0 1"} {origin ""} {weightName masses} {weightSum ""} {print 0}} {
	## applies a torque to each atom in $ids about COM	
	## currently, the weightName/weightSum options are ignored...
	variable forces
	variable coords
	variable $weightName

	## use COW (center of weight)
	switch [llength [join $origin]] {
	    3 { continue }
	    0 { foreach {origin m} [getCOW $ids $weightName] { break } }
	    default { print "WARNING: [info level 0]: origin = $origin is not a vector!; using center of weight" ;  foreach {origin m} [getCOW $ids $weighName] { break } }
	}
	set axis [vecnorm $axis];# normalize axis
	
	## build lists of information (per atom)
	set rs ""
	set rls ""
	set weights ""
	set weightSum "0"
	foreach id $ids {
	    set r [vecorthogonal [vecsub $coords($id) $origin] $axis]
	    set rl [veclength $r]
	    set w [subst $$weightName($id)]
	    
	    lappend rs $r
	    lappend rls $rl
	    lappend weights $w
	    append weightSum "+$w"
	}
	set weightSum [expr $weightSum]
	set torque 
	foreach id $ids r $rs rl $rls w $weights {
	    set force [expr $torque*$w/($rl*$weightSum)]
	    set forces($id) [vecadd $forces($id) [vecscale [vecnorm [veccross $axis $r]] $force]]
	}
    }
    proc applyTorque {ids torque {axis {0 0 1}} {origin {0 0 0}} {weightName masses} {weightSum ""} {print 0}} {
	## applies a torque to COW of $ids along $axis and about $origin
	## get center of mass
	foreach {c m} [getCOM $ids] { break }

	set axis [vecnorm $axis];# normalize axis
	set r [vecsub $c-$origin]
	set r [vecorthogonal $r $axis]
	
	## force = torque/[veclength $c-$origin]
	set force [vecnorm [veccross $axis $r]]
	set force [vecscale [expr $torque/[veclength $r]] $force ]
	
	applyForce $ids $force $weightName $weightSum $print
    }
    proc applyTorques {ids1 ids2 torque {axis "0 0 1"} {origin ""} {weightName masses} {weightSum ""} {print 0}} {
	## applies torque about the center of ids1 and ids2 so that the net force is zero
	## DOES NOT apply torque about the center of mass of each ids1/ids2
	foreach {c1 m1} [getCOM $ids1] { break }
	foreach {c2 m2} [getCOM $ids2] { break }
	
	if { [llength $origin] == 0 } { 
	    set origin [vecscale [expr {1./($m1+$m2)}] [vecadd [vecscale $m1 $c1] [vecscale $m2 $c2]]]
	}

	set axis [vecnorm $axis];# normalize axis
	set r1 [vecsub $c1 $origin]
	set r2 [vecsub $c2 $origin]
	set r1 [vecorthogonal $r1 $axis]
	set r2 [vecorthogonal $r2 $axis]
		
	## force = torque/[veclength $c-$origin]
	## want f1-f2 = 0 (already know they are pointing opposite)
	## f1 = T1/r1 ; f2 = T2/r2  ; T1 + T2 = T    =>    f1 = 
	set rl1 [veclength $r1]
	set rl2 [veclength $r2]
	
	set force [vecnorm [veccross $axis $r1]]
	set force [vecscale [expr {$torque/(${rl1}+${rl2})}] $force]
	applyForce $ids1 $force masses $m1 $print
	applyForce $ids2 $force masses -$m2 $print
    }
    proc addEnergy {name value} {
	foreach var "energyNames energy" { variable $var }
	lappend energyNames $name
	if { ![info exists energy($name)] } { 
	    set energy($name) $value
	} else { set energy($name) [expr {$energy($name)+$value}] }
    }
    proc printForce { procName name value args } {
	## prints data in a standard way; 
	variable timeStep
	variable printCount
	
	set printList "Timestep $timeStep: tclForces::$procName applied to $printCount: $name $value"
        foreach {name value} $args {
	    append printList " ; $name = $value"
	}
	print $printList
	incr printCount
    }	    
    proc printEnergy {} {
	foreach var "energyNames energy timeStep" { variable $var }
	set energyNames [lsort -unique $energyNames]
	set totalEnergy 0
	foreach name $energyNames {
	    set totalEnergy [expr {$totalEnergy + $energy($name)}]
	    print "Timestep $timeStep : $name = $energy($name) kcal/mol"
	    set energy($name) 0
	}
	if { [llength $energyNames] > 1 } {
	    print "Timestep $timeStep : totalEnergy = $totalEnergy kcal/mol"
	}
	set energyNames ""
    }

    #### internal procs (not to be exported)
    proc clearforces {} {
	foreach var "forces atoms" { variable $var }
	foreach id $atoms {
	    set forces($id) "0.0 0.0 0.0"
	}
	## also clear energies
	
    }
    proc addAtoms {} {
	variable atoms
	foreach id $atoms { ::addatom $id }
    }
    proc evalForceFiles {{ns tclForces::forces}} {
	variable forceFileContents
	namespace eval :: "namespace import ${ns}::*; $forceFileContents; namespace forget ${ns}::*"
    }			    
}

#### define ::tclForces and import things
namespace eval ::tclForces::forces {
    ## import important procs
    namespace import ::tclForces::vector::*
    namespace import ::tclForces::*

    variable exportProcs ""
}

foreach f $forceDefFiles {
    set f $tclForcesDir/$f
    if [file exists $f] {
	source $f
    } else {
	print "WARNING: force definition file $f could not be found!"
    }
}

namespace eval ::tclForces::forces {
    foreach procName [join $exportProcs] { namespace export $procName }
}


namespace eval tclForces::init {
    variable exportProcs ""
    proc append_ids { args } {
	variable ::tclForces::atoms
	lappend atoms [join $args]
    }
    # define force procedures to build a list of atoms		 
    proc printForce { procName args } {
	## prints data in a standard way;
	uplevel {
	    if { $print } {
		variable ::tclForces::printCount
		set printList "tclForces::init: \"$name\" applied (label \"$printCount\") with arguments:"
		incr printCount
	    } else {
		set printList "tclForces::init: \"$name\" applied with arguments:"
	    }
	    foreach arg $myargs {
		append printList " $arg {[subst $$arg]}"
	    }
	    print $printList
	}
    }
	
    ## create Procs that initially print what forces are being applied
    ## also add atoms!
    proc createProc { procName } {
	## this proc takes a procName = ::tclForces::forces::someProc, and creates ::tclForces::init::someProc
	set shortName [namespace tail $procName] ;# = someProc
	set oldArgList [info args $procName]
	set argList ""
	foreach arg $oldArgList {
	    if { [info default $procName $arg def] } {
		# print "setting $procName arg \"$arg {$def}\""
		set arg "$arg {$def}"
	    }
	    lappend argList "$arg"
	}
	#print "$shortName argList is $argList"
	proc $shortName $argList {
	    ## get $name of currently executing proc
	    set name [lindex [info level 0] 0]
	    set short [namespace tail $name]
	    
	    ## add atoms for all arguments of the form id.*
	    set myargs [info args $name]
	    foreach arg $myargs {
		if { [string match id* $arg] } {
		    #print "adding $arg [subst $$arg]"
		    append_ids [subst $$arg]
		}
	    }
	    printForce $name
	}
	variable exportProcs
    }
    
    ## duplicate procs in ::tclForces::forces
    #set procList [info commands ::tclForces::forces::*]
    variable exportProcs $::tclForces::forces::exportProcs
    foreach procName [join $exportProcs] {
	createProc ::tclForces::forces::$procName
	namespace export $procName
    }
    #print "init::bond args are: [info args bond]"
#     print "init::bond args are: [info args bond]"
#     foreach arg "print energy_name" {
# 	if { [info default bond $arg def] } {
# 	    print "init::bond argument $arg defaults to $def"
# 	} else { print "init::bond argument $arg has no default" }
#     }
#     print "init::restrainEachToPlane args are: [info args restrainEachToPlane]"
#     foreach arg "print energy_name" {
# 	if { [info default restrainEachToPlane $arg def] } {
# 	    print "init::restrainEachToPlane argument $arg defaults to $def"
# 	} else { print "init::restrainEachToPlane argument $arg has no default" }
#     }
}

namespace eval tclForces {
    ## do the initialization
    namespace import init::*
    
    ## parse forceFiles and store contents
    global forceFiles
    foreach file $forceFiles {
	if [file exists $file] { 
	    set ch [open $file r]
	    append forceFileContents [read $ch]
	    # 	    while {[gets $ch lines] >= 0} {
	    # 		foreach line [split $lines ";"] {
	    # 		    if {[regexp {\S*\#} $line]} { continue }
	    # 		    # set line [subst $line];# this might be bad
	    # 		    set line [namespace eval :: "subst \"$line\""]
	    # 		    eval $line
	    
	    # 		    append forceFileContents "$line ; "
	    # 		}
	    # 	    }
	    close $ch
	} else { print "WARNING: file $file does not exist" }
    }
    print "forceFileContents are: { $forceFileContents }"
    evalForceFiles tclForces::init

    ## addatoms
    set atoms [lsort -unique [join $atoms]]
    addAtoms
    clearforces
    
    ## sanity check to ensure that atoms are still present
    variable applyforce 1
    if {[llength $atoms] == 0} { set applyforce 0; print "WARNING: forces are not being applied to any atoms" } 
    
    print "Number of atoms forced: [llength $atoms]"

    ## switch force definition
    namespace forget init::*
    namespace import forces::*
    
    ## calcforces
    proc calcforces {} {
	global forcesRecalcFreq
	set vars "variable atoms masses forces coords"
	append vars " applyforce timeStep printCount forcecount"

	foreach var $vars { variable $var }
	
	if { $applyforce == 0 } { return } ;# only proceed if forces should be applied	
	# Check to see if we're about to recalculate forces
	# next time step. If so, re-add atoms
	if { $forcecount == 1 } {
	    addAtoms
	} elseif { $forcecount == 2 } {
	    ::loadcoords ::tclForces::coords
	    if { ! [array size masses] } { ::loadmasses ::tclForces::masses }
	    
	    set printCount 0
	    clearforces
	    evalForceFiles	    
	    printEnergy
	    
	    # Vital to clearconfig after calculating forces, else
	    # coordinates will be retrieved every timestep, erasing
	    # much of the potential speed gains
	    ::clearconfig

	
	} elseif { $forcecount == $forcesRecalcFreq } { set forcecount 0 }
	foreach atom $atoms {
	    addforce $atom $forces($atom)
	}
	incr timeStep
	incr forcecount
	return
    }
}
proc calcforces {} {
    # print "::calcforces time = [getstep]"
    ::tclForces::calcforces
}
