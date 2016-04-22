set targetMark    "1.00"

set targets {}
set masses {}


# Return a list with atom positions.
proc extractPdbCoords {pdbFile} {
    set r {}
    
    # Get the coordinates from the pdb file.
    set in [open $pdbFile r]
    foreach line [split [read $in] \n] {
	if {[string equal [string range $line 0 3] "ATOM"]} {
	    set x [string trim [string range $line 30 37]]
	    set y [string trim [string range $line 38 45]]
	    set z [string trim [string range $line 46 53]]
	    
	    lappend r [list $x $y $z]
	}
    }
    close $in
    return $r
}

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

set posAll [extractPdbCoords $targetAtomPdb]

# make list of atoms
set atoms {}
set pos {}
foreach target $targets {
    foreach {segname resid atom} $target { break }
    set atomindex [atomid $segname $resid $atom]
    set r [lindex $posAll [expr $atomindex - 1]]

    print "$segname $resid $atom $r"

    lappend pos $r
    lappend atoms $atomindex
    addatom $atomindex
}

set numatoms [llength $atoms]

if { $numatoms <=0 } {
    print "WARNING: no target atoms have been detected"
}

print "Number of atoms forced: $numatoms"
print "positionZ0: $positionZ0"
print "spring: $spring"
print "restK: $restK"
print "lz: $lz"

# Zero out initial forces
set forces {}
foreach index $atoms {
    lappend forces "0.0 0.0 0.0"
}

set forcecount $forcesRecalcFreq
set printcount 0
set kcalPerMolA 69.4769460332
set halfZ [expr 0.5*$lz]

proc calcforces { } {
    global atoms numatoms masses forces
    global forcesRecalcFreq
    global spring positionZ0
    global forcecount printcount
    global kcalPerMolA pos restK halfZ lz
    
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
	
	# Handle periodic boundary conditions.
	set dz [expr $positionZ0 - $z]
	if {$dz < -$halfZ} {
	    set dz [expr $dz + $lz]
	} elseif {$dz >= $halfZ} {
	    set dz [expr $dz - $lz]
	}
	set fz [expr $spring*$dz]
	
	# Now calculate forces
	set forces {}
	foreach atom $atoms p $pos {
	    set fx [expr $restK*([lindex $p 0] - [lindex $coords($atom) 0])]
	    set fy [expr $restK*([lindex $p 1] - [lindex $coords($atom) 1])]
	    lappend forces [list $fx $fy $fz]
	}
	
	set ft [expr $kcalPerMolA*$numatoms*$fz]
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
