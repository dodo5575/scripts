# Calculate the center of mass for a trajectory.
# Writes step time(ns) length (A) angle(rad) to a file.
# to use: vmd -dispdev text -e trackPositionDcd.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set dcdFreq 5000
set selText "segname HAIR"
set timestep 1.0
set stride 20
set resIdList {127 126 125 124 123 122 121 120 119 118 103 104 105 106 107 108 109 110 111}
#set resIdList {125 104}

# Input:
set psf pore+dna-all.psf
set pdb pore+dna-all.pdb
set dcdPrefix /Scr/nanopore/jcomer/myhairpin/loop_first_dcd/loop_first_2V
set dcdSuffix ".dcd"
set dcdSet {0 1 2 3}
# Output:
set outPrefix backbone_loop_first_2V

# Get the time change between frames in nanoseconds.
set dt [expr 1.0e-6*$timestep*$dcdFreq*$stride]
set dcdEnd [llength $dcdSet]

# Open the output files.
set out {}
foreach res $resIdList {
    lappend out [open ${outPrefix}_${res}.dat w]
}

# Load the system.
mol load psf $psf pdb $pdb

# Get the atom indices for atoms over which the distance will be measured.
set bbName {P O5' C5' C4' C3' O3' P}
set bbRes {0 0 0 0 0 0 1}
set bbIndex0 {}
set bbIndex1 {}
foreach res $resIdList {
    set ind0 {}
    set ind1 {}
    
    for {set i 1} {$i < [llength $bbName]} {incr i} {
	set r0 [lindex $bbRes [expr $i-1]]
	set n0 [lindex $bbName [expr $i-1]]
	set r1 [lindex $bbRes $i]
	set n1 [lindex $bbName $i]
	
	set sel0 [atomselect top "($selText) and resid [expr $res+$r0] and name $n0"]
	set sel1 [atomselect top "($selText) and resid [expr $res+$r1] and name $n1"]
	lappend ind0 [lindex [$sel0 get index] 0]
	lappend ind1 [lindex [$sel1 get index] 0]
	
	$sel0 delete
	$sel1 delete
    }
    
    lappend bbIndex0 $ind0
    lappend bbIndex1 $ind1
}
puts $bbIndex0
puts $bbIndex1

# Loop over the dcd files.
set nFrames0 1
for {set dcd 0} {$dcd < $dcdEnd} {incr dcd} {
    # Load the trajectory.
    animate delete all
    set dcdNum [lindex $dcdSet $dcd]
    mol addfile "${dcdPrefix}${dcdNum}${dcdSuffix}" type dcd step $stride waitfor all
    set nFrames [molinfo top get numframes]
    puts [format "Reading %i frames." $nFrames]

    # Start at "startFrame" and move forward, computing
    # the sasa at each step.
    for {set f 0} {$f < $nFrames} {incr f} {
	molinfo top set frame $f
	
	# Get the time in nanoseconds for this frame.
	#set step [expr ($nFrames0+$f)*$dcdFreq*$stride]
	set t [expr ($nFrames0+$f)*$dt]

	set outText {}
	lappend $outText $t

	foreach res $resIdList o $out ind0 $bbIndex0 ind1 $bbIndex1 {
	    # Determine the length of the backbone.
	    set length 0.0
	    foreach i0 $ind0 i1 $ind1 {
		set length [expr $length + [measure bond [list $i0 $i1]]]
	    }

	    # Write the length and angle.
	    puts $o "$t $length"
	}
	
	puts -nonewline [format "FRAME %i: " $f]
	puts "$t $length"
    }
    set nFrames0 [expr $nFrames+$nFrames0]
}

foreach o $out {
    close $o
}

exit



