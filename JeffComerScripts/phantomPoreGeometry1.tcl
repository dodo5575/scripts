# phantomPoreMembrane.tcl
# Force a selection out of double cone pore in a membrane.
# Inputs: poreForce forcesRecalcFreq targetAtomPdb poreGeometryFile
# The largest radii must be on the top and bottom!
print "Starting tclForces: phantomPoreGeometry"

# Define the target beta value.
set targetMark "1.00"

# Read the pore radius as a function of z from the file.
proc readGeometry {fileName} {
    set z {}
    set s {}

    set in [open $fileName r]
    foreach line [split [read $in] \n] {
	set tok [concat $line]
	
	if {[llength $tok] == 2} {
	    lappend z [lindex $tok 0]
	    lappend s [lindex $tok 1]
	}
    }
    close $in

    return [list $z $s]
}

# Find the atoms to which tclForces will be applied by
# checking their beta column.
proc readTargetAtoms {fileName targetMark} {
    set atomList {}

    set in [open $fileName r]
    foreach line [split [read $in] \n] {
	set record [string trim [string range $line 0 5]]
	set serial [string trim [string range $line 6 10]]
	set name [string trim [string range $line 12 15]]
	set altLoc [string trim [string range $line 16 16]]
	set resName [string trim [string range $line 17 19]]
	set chainId [string trim [string range $line 21 21]]
	set resId [string trim [string range $line 22 25]]
	set iCode [string trim [string range $line 26 26]]
	set x [string trim [string range $line 30 37]]
	set y [string trim [string range $line 38 45]]
	set z [string trim [string range $line 46 53]]
	set occupancy [string trim [string range $line 54 59]]
	set beta [string trim [string range $line 60 65]]
	set segName [string trim [string range $line 72 75]]
	set element [string trim [string range $line 76 77]]
	set charge [string trim [string range $line 78 79]]

	if {[string match "ATOM*" $record] || [string match "HETA*" $record]} {
	    if {[string trim $beta] == $targetMark} {
		lappend atomList [list $segName $resId $name]
	    }
	}
    }
    close $in
           
    return $atomList
}

# Extract the geometry of the pore.
set geo [readGeometry $poreGeometryFile]
set poreZ [lindex $geo 0]
set poreS [lindex $geo 1]
set sFirst [lindex $poreS 0]
set sLast [lindex $poreS end]

# Extract the z interval.
set dz [expr [lindex $poreZ 1]-[lindex $poreZ 0]]
set z0 [lindex $poreZ 0]
set z1 [lindex $poreZ end]
set zMid [expr 0.5*($z0+$z1)]

# Make a list of pore slopes and normal forces.
set normalS {}
set normalZ {}
set slope {}
set s1 [lindex $poreS 0]
for {set i 1} {$i < [llength $poreS]} {incr i} {
    set s0 $s1
    set s1 [lindex $poreS $i]

    lappend slope [expr ($s1-$s0)/$dz]
    set nz [expr $s1-$s0]
    set ns [expr -$dz]
    set mag [expr sqrt($nz*$nz + $ns*$ns)] 
    lappend poreForceZ [expr $poreForce*$nz/$mag]
    lappend poreForceS [expr $poreForce*$ns/$mag]
}

# Make a list of indices.
set atomList [readTargetAtoms $targetAtomPdb $targetMark]
set targetAtoms {}
foreach atom $atomList {
    foreach {segName resId name} $atom {break}
    set index [atomid $segName $resId $name]
    addatom $index
    lappend targetAtoms $index 
    #print "TCLFORCESATOM $index $segName $resId $name"
}
if {[llength $targetAtoms] == 0} {
	print "WARNING: No atoms to which to apply tclForces."
}
print "Registered [llength $targetAtoms] target atoms."

# Returns the force if the atom is within the pore wall,
# a scalar zero otherwise.
proc poreForce {x y z} {
    global poreForce slope poreForceS poreForceZ
    global poreS poreZ sFirst sLast z0 dz zMid
    
    # Find the nodes that immediately surround this point.
    set j0 [expr int(floor(($z-$z0)/$dz))]
    set j1 [expr $j0 + 1]
        
    # Quit if we're above or below the pore.
    if {[expr $j0 < 0 || $j1 >= [llength $poreZ]]} {
	return 0
    }

    # Find the radii.
    set s [expr sqrt($x*$x + $y*$y)]
        
    if {[expr $s < [lindex $poreS $j0]]} {
	return 0
    }

    # Find the radius of the pore at this position.
    set s0 [lindex $poreS $j0]
    set m [lindex $slope $j0]
    set sp [expr $s0 + $m*($z-$dz*$j0)]
    # Is the atom away from the pore walls?
    if {[expr $s < $sp]} {
	return 0
    }
    
    # Exert a force from the pore walls. 
    set fs [lindex $poreForceS $j0]
    return [list [expr $fs*$x/$s] [expr $fs*$y/$s] [lindex $poreForceZ $j0]]
}

###################################################################
# This procedure is executed at each time step.
set forceAtoms {}
set forces {}
set forcesRecalcFreq1 [expr $forcesRecalcFreq-1]
set recalcCount forcesRecalcFreq
print "Starting calcforces..."
print "Force sums: step number_of_atoms force_sum"
print "Forces: step index x y z Fx Fy Fz"
proc calcforces {} {
    global recalcCount forcesRecalcFreq forcesRecalcFreq1
    global targetAtoms
    global forceAtoms forces
		
    if {$recalcCount == $forcesRecalcFreq1} {
	# Get the atom coordinates for the next force calculation.
	foreach index $targetAtoms {
	    addatom $index
	}
    } elseif {$recalcCount >= $forcesRecalcFreq} {
	# Recalculate the forces.
	set forceSum [list 0.0 0.0 0.0]
	loadcoords coords
	set step [getstep]
				
	set forceAtoms {}
	set forces {}
	foreach index $targetAtoms {
	    foreach {x y z} $coords($index) {break}
	    # Compute the force, if there is one.
	    set f [poreForce $x $y $z]
			
	    if {[llength $f] != 1} {
		lappend forceAtoms $index
		lappend forces $f
		#set forceSum [vecadd $forceSum $f]
	    }
	}
	
	print "TCLFORCESSUM $step [llength $forces]"
	clearconfig
	set recalcCount 0
    }
	
    # Add the forces.
    foreach i $forceAtoms f $forces {
	addforce $i $f
    }
    
    incr recalcCount
    return
}



