# Calculate the center of mass for a trajectory.
# Writes step time(ns) and position (nm) to a file.
# to use: vmd -dispdev text -e trackPositionDcd.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set dcdFreq 5000
set selText "ions and abs(z) < 96"
set timestep 1.0
set stride 1

# Input:
set psf nw_pore+dna-all.psf
set pdb nw_pore+dna-all.pdb
set dcdPrefix /Scr/nanopore/jcomer/myhairpin/conform4_dcd/nw_
set dcdSuffix "_4V.dcd"
set dcdSet {run1 run2 run3 run4 run5 run6 fix0 fix1 fix2 fix3 fix4 fix5 fix6 fix7}
# Output:
set outPrefix ion_conform4_4V[lindex $dcdSet 0]-[lindex $dcdSet end].dat

# Get the time change between frames in nanoseconds.
set dt [expr 1.0e-6*$timestep*$dcdFreq*$stride]
set dcdEnd [llength $dcdSet]

# Open the output file.
set out [open $outPrefix w]

# Load the system.
mol load psf $psf pdb $pdb

# Loop over the dcd files.
set nFrames0 0
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

	# Get the number in the selection.
	set sel [atomselect top $selText]
	set n [$sel num]
	$sel delete

	# Write the time and position.
	puts $out "$t $n"
	
	puts -nonewline [format "FRAME %i: " $f]
	puts "$t $n"
    }
    set nFrames0 [expr $nFrames+$nFrames0]
}

close $out
exit



