# Compute the solvent-accessible surface area for a trajectory.
# to use: vmd -dispdev text -e sasa.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

# SASA radius
set radius 2.0
set startFrame 0
set stride 1
set selText "protein"
set dcdFreq 1000

# Input:
set pdb  1p2y.pdb
set psf  1p2y.psf
set dcd  run10_100V.dcd
# Output:
set outFile water_heme_run10.txt

# Make the output in picoseconds.
set rate [expr $dcdFreq/1000.]

# Load the system.
mol load psf $psf pdb $pdb
set sel [atomselect top $selText]

# Load the trajectory.
animate delete all
mol addfile $dcd waitfor all
set nFrames [molinfo top get numframes]
puts [format "Reading %i frames." $nFrames]

# Open the output file.
set out [open $outFile w]
puts $out "SASA for $psf with trajectory $dcd and radius $radius"
puts $out "t(ps) SASA(A^2)"

# Start at "startFrame" and move forward "stride" frames, computing
# the SASA at each step.
set n 1
for {set f $startFrame} {$f < $nFrames && $n > 0} {incr f $stride} {
	molinfo top set frame $f
	puts [format "FRAME %i" [molinfo top get frame]]

	set sasa [measure sasa $radius $sel]
	set t [expr $f*$rate]
	puts $out "$t $sasa"
}
close $out
exit




