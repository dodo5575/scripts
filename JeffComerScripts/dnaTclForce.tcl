set targetMark    "1.00"

set targets {}
set masses {}

set inStream [open $targetAtomPdb r]
foreach line [split [read $inStream] \n] {
    set string1 [string range $line 0 3]
    set string2 [string range $line 6 10]
    set string3 [string range $line 17 20]	
    set string4 [string range $line 12 15]
    set string5 [string range $line 46 53]
    set string6 [string range $line 72 75]
    set string673 [string range $line 73 75]
    set string7 [string range $line 22 25]
    set string8 [string range $line 62 65]
    set string9 [string range $line 55 60]
    
    
    if { ([string equal $string1 {ATOM}] || \
 	      [string equal $string1 {HETA}] ) && \
 	     [string equal $targetMark $string8] } {	
 	lappend targets "[string trim $string6]\
 			    [string trim $string7] [string trim $string4]"
	lappend masses "[string trim $string9]"
    }
}
close $inStream


# make list of atoms
set atoms {}
foreach target $targets {
    foreach {segname resid atom} $target { break }
    print "$segname $resid $atom"
    set atomindex [atomid $segname $resid $atom]
    lappend atoms $atomindex
    addatom $atomindex
}

set numatoms [llength $atoms]

if { $numatoms > 0 } {
    set applyforce 1
} else {
    print "WARNING: no target atoms have been detected"
    set applyforce 0
}

# Take force factor from NAMD config file
set forcemult [expr double($forceK)/$numatoms]
    
print "Number of atoms forced: $numatoms"
print "Force multiplier (per atom): $forcemult"

# Zero out initial forces
set forces {}
foreach index $atoms {
    lappend forces "0.0 0.0 0.0"
}

set forcecount $forcesRecalcFreq
set printcount 0
set kcalPerMolA 69.4769460332

proc calcforces { } {
    global atoms numatoms forcemult masses  forces
    global forcesRecalcFreq
    global spring positionZ0
    global forcecount printcount
    global kcalPerMolA
    
    # Apply forces
    foreach atom $atoms force $forces {
	addforce $atom $force
    }
    
    # Check to see if we're about to recalculate forces
    # next time step. If so, clearconfig and re-add atoms
    if { $forcecount == [expr $forcesRecalcFreq - 1] } {
	#    print "Adding atoms prior to reconfiguring forces at $printcount"
	foreach atom $atoms {
	    addatom $atom
	}
    }
    
    if { $forcecount == $forcesRecalcFreq } {
	#    print "Recalculating forces at $printcount"
	
	# Get coordinates
	loadcoords coords
	
	# First calculate center
	set coordsum "0.0 0.0 0.0"
	foreach atom $atoms {
	    set coordsum [vecadd $coordsum $coords($atom)]
	}
	set center [vecscale [expr 1.0/$numatoms] $coordsum]
	#print "Center = $center"
	set x [lindex $center 0]
	set y [lindex $center 1]
	set z [lindex $center 2]
	set force "[expr -$forcemult*$x] [expr -$forcemult*$y] [expr -$spring*($z-$positionZ0)]"
	
	# Now calculate forces
	set forces {}
	foreach atom $atoms mass $masses {
	    lappend forces $force
	}
	
	set ft [expr $kcalPerMolA*$numatoms*[lindex $force 2]]
	print "TCLDNAFORCE [getstep] $center $ft"
	
	set forcecount 0
	
	# Vital to clearconfig after calculating forces, else
	# coordinates will be retrieved every timestep, erasing
	# much of the potential speed gains
	clearconfig
    }
    incr forcecount

    incr printcount
    return
}
