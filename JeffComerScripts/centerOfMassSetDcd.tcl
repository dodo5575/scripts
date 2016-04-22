# Calculate the center of mass for a trajectory.
# Writes step time(ns) and position (nm) to a file.
# to use: vmd -dispdev text -e trackPositionDcd.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set dcdFreq 5000
set resSet {25 26 27 28 29 30}
set selText0 "segname ADNA"
set selName0 ""
set timestep 1.0
set stride 1

# Input:
set voltage 8
set psf pore+dna-all.psf
set pdb pore+dna-all.pdb
set dcdPrefix /Scr/nanopore/jcomer/double/trap2.0_dcd/trap2.0_pulse200ps_${voltage}V
set dcdSuffix ".dcd"
set dcdSet {2}
# Output:r
set outDir analysis
set outName p2.0_pulse${voltage}V[lindex $dcdSet end]

# Generate the selections.
set selText {}
set selName {}
foreach res $resSet {
    lappend selName "${selName0}${res}"
    lappend selText "$selText0 and resid $res"
}

# Get the time change between frames in nanoseconds.
set dt [expr 1.0e-6*$timestep*$dcdFreq*$stride]
set dcdEnd [llength $dcdSet]

# Open the output file.
set out {}
foreach sn $selName {
    lappend out [open $outDir/comz_${sn}_${outName}.dat w]
}

# Load the system.
mol load psf $psf pdb $pdb
set sel {}
foreach st $selText {
    lappend sel [atomselect top $st]
}

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
	puts -nonewline [format "FRAME %i: " $f]
	puts -nonewline "$t"

	foreach s $sel o $out {
	    set r [measure center $s weight mass]
	    # Get the center of mass position in nm for this frame.
	    set z [expr [lindex $r 2]/10.0]

	    # Write the time and position.
	    puts $o "$t $z"
	
	    puts -nonewline " $z"
	}
	puts ""
    }
    set nFrames0 [expr $nFrames+$nFrames0]
}

foreach o $out s $sel {
    close $o
    $s delete
}
exit



