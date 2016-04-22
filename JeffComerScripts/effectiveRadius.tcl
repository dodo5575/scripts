# Calculate the center of mass for a trajectory.
# Writes time(ns) position(nm) to a file.
# to use: vmd -dispdev text -e trackPositionDcd.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set dcdFreq 5000
set selText "segname HAIR"
set stride 1
set zRange 10.0

# Input:
set pdb  pore+dna-all.pdb
set psf  pore+dna-all.psf
set dcdPrefix  /Projects/alek/HAIRPIN/20_HAIR-PORE/5_RUNS/hairpin-E6V-
set dcdSuffix ".dcd"
set dcdSet {1 2 3 4 5}
# Output:
set outFile erad_pore2.0.txt

# Get the time change between frames in femtoseconds.
set dcdEnd [llength $dcdSet]

# Open the output file.
set out [open $outFile w]

# Load the system.
mol load psf $psf pdb $pdb
set sel [atomselect top $selText]

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

	set position [$sel get {x y z}]
	
	foreach r $position {
		foreach {x y z} $r {break}
		set s [expr sqrt($x*$x + $y*$y)]
		if {[expr abs($z) < $zRange]} {
			puts $out "$z $s"
		}
	}
	
	puts "frame $f"
}
set nFrames0 [expr $nFrames+$nFrames0]
}

close $out
$sel delete
exit




