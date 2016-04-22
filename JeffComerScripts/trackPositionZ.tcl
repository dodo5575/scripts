# Calculate position of the center of mass of the selection
# for a trajectory.
# to use: vmd -dispdev text -e trackPositionZ.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set dcdFreq 250000
set selText "segname HAIR"
set startFrame 0
set timestep 1.0

# Input:
set pdb  translocate.pdb
set psf  translocate.psf
set dcd  translocate.dcd
# Output:
set outFile pos_6V.dat

# Get the time change between frames in femtoseconds.
set dt [expr $timestep*$dcdFreq]

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

# Start at "startFrame" and move forward, computing
# current at each step.
set n 1
for {set f $startFrame} {$f < $nFrames && $n > 0} {incr f} {
	molinfo top set frame $f
	
	set r [measure center $sel weight mass]
	set z [lindex $r 2]
	
	# Get the time in nanoseconds for this frame.
	set t [expr ($f+0.5)*$dt*1.e-6]
	
	puts $out "$t $z"
	puts -nonewline [format "FRAME %i: " $f]
	puts "$t $z"
}
close $out
mol delete top
exit



