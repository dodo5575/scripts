# Calculate the center of mass for a trajectory.
# Writes time(ns) position(nm) to a file.
# to use: vmd -dispdev text -e trackPositionDcd.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set dcdFreq 5000
set selText "resname TIP3 and name OH2"
set stride 1
set zRange 1.0

# Input:
set pdb  pore1_6_all.pdb
set psf  pore1_6_all.psf
set dcdPrefix  p1.6_open
set dcdSuffix ".dcd"
set dcdSet {1 2 3 4}
# Output:
set outFile wrad_pore1.6.txt

# Get the time change between frames in femtoseconds.
set dcdEnd [llength $dcdSet]

# Open the output file.
set out [open $outFile w]

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

	set sel [atomselect top "abs(z) < $zRange and ($selText)"]
	set position [$sel get {x y z}]
	$sel delete
	
	foreach r $position {
		foreach {x y z} $r {break}
		set s [expr sqrt($x*$x + $y*$y)]
		puts $out "$z $s"
	}
	
	puts "frame $f"
}
set nFrames0 [expr $nFrames+$nFrames0]
}

close $out
exit




