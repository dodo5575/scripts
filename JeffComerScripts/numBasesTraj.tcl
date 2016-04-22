# Calculate the distance between two CoMs for a trajectory.
# Writes step time(ns) position (nm) cross-sectional radius (nm) to a file.
# to use: vmd -dispdev text -e stretchDcd.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

# Requires: stride timestep dcdFreq psf pdb xsc dcd outSuffix moiety displayPeriod
switch $moiety {
    hpDNA {set selText "segname HAIR"}
    dsDNA {set selText "segname ADNA BDNA"}
    default {
	puts stderr "ERROR: Unrecognized moiety $moiety."
	exit
    }
}

# Get the time change between frames in nanoseconds.
set dt [expr 1.0e-6*$timestep*$dcdFreq*$stride]
set dcdEnd [llength $dcd]

# Open the output file.
set out [open output/bases${dz}_${outSuffix} w]

# Load the system.
mol load psf $psf pdb $pdb

# Loop over the dcd files.
set nFrames0 0
foreach dcdFile $dcd {
    # Load the trajectory.
    animate delete all
    mol addfile $dcdFile type dcd step $stride waitfor all
    set nFrames [molinfo top get numframes]
    puts [format "Reading %i frames." $nFrames]

    # Move forward computing the center-of-mass at every step.
    for {set f 0} {$f < $nFrames} {incr f} {
	molinfo top set frame $f

	# Get the time in nanoseconds for this frame.
	set t [expr ($nFrames0+$f)*$dt]
	
	set sel [atomselect top "($selText) and abs(z) < $dz"]
	set n [llength [lsort -unique [$sel get {segname resid}]]]
	$sel delete

	# Write the time and number.
	puts $out "$t $n"
		
	if {$f % $displayPeriod == 0} {
	    puts -nonewline [format "FRAME %i: " $f]
	    puts "$t $n"
	}
    }
    set nFrames0 [expr $nFrames+$nFrames0]
}
mol delete top
close $out



