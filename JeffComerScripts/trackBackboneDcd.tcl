# Calculate the center of mass for a trajectory.
# Writes step time(ns) length (A) angle(rad) to a file.
# to use: vmd -dispdev text -e trackPositionDcd.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set dcdFreq 5000
set selText "segname HAIR"
set timestep 1.0
set stride 1
set cutoff 3.0
set resIdList {127 126 125 124 123 122 121 120 119 118 103 104 105 106 107 108 109 110 111}
set up {0.0 0.0 1.0}

# Input:
set pdb tall_sys1.pdb
set psf tall_sys1.psf
set dcdPrefix uz1_zip
set dcdSuffix ".dcd"
set dcdSet {0 1 2 3}
# Output:
set outPrefix backbone_uz1

# Get the time change between frames in femtoseconds.
set dt [expr $timestep*$dcdFreq]
set dcdEnd [llength $dcdSet]

# Open the output file.
set out {}
foreach res $resIdList {
    lappend out [open ${outPrefix}_${res}.dat w]
}

# Load the system.
mol load psf $psf pdb $pdb

# Loop over the dcd files.
set nFrames0 1
for {set dcd 0} {$dcd < $dcdEnd} {incr dcd} {
    # Load the trajectory.
    animate delete all
    set dcdNum [lindex $dcdSet $dcd]
    mol addfile "${dcdPrefix}${dcdNum}${dcdSuffix}" waitfor all
    set nFrames [molinfo top get numframes]
    puts [format "Reading %i frames." $nFrames]

    # Start at "startFrame" and move forward, computing
    # the sasa at each step.
    for {set f 0} {$f < $nFrames} {incr f $stride} {
	molinfo top set frame $f
	
	# Get the time in nanoseconds for this frame.
	set step [expr ($nFrames0+$f)*$dcdFreq]
	set t [expr ($nFrames0+$f)*$dt*1.e-6]

	set outText {}
	lappend $outText $t
	set angle0 0.0
	set angle1 0.0

	foreach res $resIdList o $out {
	    set text0 "($selText) and resid [expr $res-1] O5'"
	    set text "($selText) and resid $res and name O5'"
	    
	    # Determine the length and angle of the backbone.
	    set sel0 [atomselect top $text0]
	    set sel [atomselect top $text]
	    
	    if {[$sel0 num] == 0 || [$sel num] == 0} {
		puts "Cannot select backbone for residue $res."
		exit
	    }
	    
	    set r0 [lindex [$sel0 get {x y z}] 0]
	    set r [lindex [$sel get {x y z}] 0]
	    set d [vecsub $r $r0]
	    set length [veclength $d]
	    set comp [vecscale [expr 1./$length] [vecdot $d $up]]
	    set angle [expr acos($comp)]
	    $sel0 delete
	    $sel delete

	    # Write the length and angle.
	    puts $o "$step $t $length $angle"
	}
	
	puts -nonewline [format "FRAME %i: " $f]
	puts "$step $t $length $angle"
    }
    set nFrames0 [expr $nFrames+$nFrames0]
}

foreach o $out {
    close $o
}

exit



