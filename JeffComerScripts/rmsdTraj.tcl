# Calculate the current for a trajectory.
# Results are in "time(ns) current(nA)"
# to use: vmd -dispdev text -e electricCurrentZ.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set dcdFreq 16000
set selText "protein"
set startFrame 0
set timestep 1.0

# Input:
set dir supplementary-files/thermal-stability
set prefix bpti-reduced-310K
set psf ${dir}/bpti-reduced-system.psf
set dcd ${dir}/${prefix}.dcd
# Output:
set outFile rmsd_${prefix}.dat

# Get the time change between frames in femtoseconds.
set dt [expr $timestep*$dcdFreq]

# Load the reference.
set refMol [mol load psf $psf]
set nFrames [molinfo top get numframes]
puts [format "Reading %i frames." $nFrames]
mol addfile $dcd waitfor all first 0 last 0
set ref [atomselect top $selText]

# Load the system.
mol load psf $psf
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
for {set f [expr $startFrame]} {$f < $nFrames} {incr f} {
	molinfo top set frame $f
		
	set currentZ [measure rmsd $ref $sel]
	# Get the time in nanoseconds for this frame.
	set t [expr ($f+1)*$dt*1.e-6]
			
	# Write the current.
	puts $out "$t $currentZ"
	puts -nonewline [format "FRAME %i: " $f]
	puts "$t $currentZ"
}
close $out
$ref delete
$sel delete
mol delete top
mol delete top
exit



