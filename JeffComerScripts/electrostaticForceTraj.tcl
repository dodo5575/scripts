# Calculate the distance between two CoMs for a trajectory.
# Writes step time(ns) position (nm) cross-sectional radius (nm) to a file.
# to use: vmd -dispdev text -e stretchDcd.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

# Requires: stride timestep dcdFreq psf pdb xsc dcd outSuffix moiety displayPeriod
set membrane 99
switch $moiety {
    hpDNA {
	set selText "segname HAIR"
    }
    dsDNA {
	set selText "segname ADNA BDNA"
    }
    default {
	puts stderr "ERROR: Unrecognized moiety $moiety."
	exit
    }
}

# Get the time change between frames in nanoseconds.
set dt [expr 1.0e-6*$timestep*$dcdFreq*$stride]
set dcdEnd [llength $dcd]

# Open the output files.
set out [open output/fz_high_${outSuffix} w]

# Load the system.
mol load psf $psf pdb $pdb
mol addFile $grid waitfor all
set sel [atomselect top $selText]
foreach null {0} {set charge [$sel get charge]}

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

	# Get the force in piconewtons for this frame.
	set fz 0.0
	set eFieldZ [$sel get interpvol0]
	set pos [$sel get z]
	foreach ez $eFieldZ q $charge z $pos {
	    if {abs($z) > $membrane} {continue}
       	    if {[string equal nan $ez]} {set ez 0.0}
	    set fz [expr $fz + $ez*$q]
	}
	set fz [expr 1602.176487*$fz]

	# Write the time (ns) and force (pN).
	puts $out "$t $fz"
	
	if {$f % $displayPeriod == 0} {
	    puts -nonewline [format "FRAME %i: " $f]
	    puts "$t $fz"
	}
    }
    set nFrames0 [expr $nFrames+$nFrames0]
}
$sel delete
mol delete top

close $out



