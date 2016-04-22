# Calculate the center of mass for a trajectory.
# Writes step time(ns) and position (nm) to a file.
# to use: vmd -dispdev text -e trackPositionDcd.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set dcdFreq 5000
set selText "segname HAIR and resid 114 115"
set timestep 1.0
set stride 20

# Input:
set psf pore+dna-all.psf
set pdb pore+dna-all.pdb
set dcdPrefix /Scr/nanopore/jcomer/myhairpin/loop_first_dcd/hairpin_first_1V
set dcdSuffix ".dcd"
set dcdSet {7 8 9}
# Output:
set outPrefix hairpin_first_1V[lindex $dcdSet 0]-[lindex $dcdSet end].txt

# Get the time change between frames in nanoseconds.
set dt [expr 1.0e-6*$timestep*$dcdFreq*$stride]
set dcdEnd [llength $dcdSet]

# Open the output file.
set out [open $outPrefix w]

# Load the system.
mol load psf $psf pdb $pdb
set sel [atomselect top $selText]

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

	set r [measure center $sel weight mass]
	# Get the center of mass position in nm for this frame.
	set z [expr [lindex $r 2]/10.0]

	# Write the time and position.
	puts $out "$t $z"
	
	puts -nonewline [format "FRAME %i: " $f]
	puts "$t $z"
    }
    set nFrames0 [expr $nFrames+$nFrames0]
}

close $out
exit



