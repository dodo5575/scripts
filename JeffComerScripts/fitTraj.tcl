# Calculate the distance between two CoMs for a trajectory.
# Writes step time(ns) position (nm) cross-sectional radius (nm) to a file.
# to use: vmd -dispdev text -e stretchDcd.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

# Requires: stride timestep dcdFreq psf pdb xsc dcd outSuffix moiety displayPeriod
switch $moiety {
    hpDNA {
	set selText "segname HAIR and resid 102 to 111 118 to 127 and backbonetype nucleicback"
	set psf0 wad_all.psf
	set pdb0 helix_0V0.pdb
    }
    dsDNA {
	
	set selText "(segname ADNA and resid 1 to 10) or (segname BDNA and resid 102 to 111) and backbonetype nucleicback"
	set psf0 pore+dna-all.psf
	set pdb0 pore+dna-all.pdb
    }
    default {
	puts stderr "ERROR: Unrecognized moiety $moiety."
	exit
    }
}

# Get the time change between frames in nanoseconds.
set dt [expr 1.0e-6*$timestep*$dcdFreq*$stride]
set dcdEnd [llength $dcd]

# Open the output file.
set out [open output/fit_${outSuffix} w]

# Load the comparison system.
mol load psf $psf0 pdb $pdb0
set sel0 [atomselect top $selText]

# Load the system.
mol load psf $psf pdb $pdb
set sel [atomselect top $selText]

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
	
	# Fit and find the minimum RMSD (nm).
	set fit [measure fit $sel $sel0]
	$sel move $fit
	set d [expr 0.1*[measure rmsd $sel $sel0]]

	# Write the time and rmsd.
	puts $out "$t $d"
		
	if {$f % $displayPeriod == 0} {
	    puts -nonewline [format "FRAME %i: " $f]
	    puts "$t $d"
	}
    }
    set nFrames0 [expr $nFrames+$nFrames0]
}

$sel0 delete
$sel delete
mol delete top
mol delete top
close $out



