# Calculate the distance between two CoMs for a trajectory.
# Writes step time(ns) and position (nm) to a file.
# to use: vmd -dispdev text -e stretchDcd.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set dcdFreq 5000
set selText0 "segname HAIR and resid 102 127"
set selText1 "segname HAIR and resid 111 118"
set timestep 1.0
set stride 20

# Input:
set psf pore+dna-all.psf
set pdb pore+dna-all.pdb
set dcdPrefix /Scr/nanopore/jcomer/myhairpin/loop_first_dcd/loop_first_1.5V
set dcdSuffix ".dcd"
set dcdSet {0 1 2 3 4 5 6}
# Output:
set outPrefix str_loop_first_1.5V[lindex $dcdSet 0]-[lindex $dcdSet end].dat

# Get the time change between frames in nanoseconds.
set dt [expr 1.0e-6*$timestep*$dcdFreq*$stride]
set dcdEnd [llength $dcdSet]

# Open the output file.
set out [open $outPrefix w]

# Load the system.
mol load psf $psf pdb $pdb
set sel0 [atomselect top $selText0]
set sel1 [atomselect top $selText1]

# Loop over the dcd files.
set nFrames0 1
for {set dcd 0} {$dcd < $dcdEnd} {incr dcd} {
    # Load the trajectory.
    animate delete all
    set dcdNum [lindex $dcdSet $dcd]
    mol addfile "${dcdPrefix}${dcdNum}${dcdSuffix}" type dcd step $stride waitfor all
    set nFrames [molinfo top get numframes]
    puts [format "Reading %i frames." $nFrames]

    # Move forward computing the center-of-mass at every step.
    for {set f 0} {$f < $nFrames} {incr f} {
	molinfo top set frame $f

	# Get the time in nanoseconds for this frame.
	set t [expr ($nFrames0+$f)*$dt]

	set r0 [measure center $sel0 weight mass]
	set r1 [measure center $sel1 weight mass]
	# Get the distance in nanometers.
	set d [expr 0.1*[veclength [vecsub $r1 $r0]]]

	# Write the time and position.
	puts $out "$t $d"
	
	puts -nonewline [format "FRAME %i: " $f]
	puts "$t $d"
    }
    set nFrames0 [expr $nFrames+$nFrames0]
}

close $out
exit



