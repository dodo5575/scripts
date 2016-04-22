# Calculate dipole moment of the selection
# for a trajectory.
# to use: vmd -dispdev text -e dipoleMomentZDiff.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set constraint 1.0
set dcdFreq 100
set selText "all"
set startFrame 4
set timestep 1.0

# Input:
set psf  ../1_build/membrane_hex.psf
set dcd field${constraint}.dcd
set dcd0 null${constraint}.dcd
# Output:
set outFile dipole${constraint}.dat

# Get the time change between frames in femtoseconds.
set dt [expr $timestep*$dcdFreq]

# Load the system.
set traj [mol load psf $psf dcd $dcd]
set sel [atomselect $traj $selText]
set traj0 [mol load psf $psf dcd $dcd0]
set sel0 [atomselect $traj0 $selText]

# Choose nFrames to be the smaller of the two.
set nFrames [molinfo $traj get numframes]
set nFrames0 [molinfo $traj0 get numframes]
if {$nFrames0 < $nFrames} {set nFrames $nFrames}
puts [format "Reading %i frames." $nFrames]

# Open the output file.
set out [open $outFile w]

# Start at "startFrame" and move forward, computing
# the dipole moment at each step.
set sum 0.
set sumSq 0.
set n 1
puts "t (ns)\tp_z (e A)\tp0_z (e A)\tp_z-p0_z (e A)"
for {set f $startFrame} {$f < $nFrames && $n > 0} {incr f} {
	$sel frame $f
	$sel0 frame $f
	
	# Obtain the dipole moment along z.
	set p [measure dipole $sel]
	set p0 [measure dipole $sel0]
	set z [expr [lindex $p 2] - [lindex $p0 2]]
	
	# Get the time in nanoseconds for this frame.
	set t [expr ($f+0.5)*$dt*1.e-6]
	
	puts $out "$t $z"
	puts -nonewline [format "FRAME %i: " $f]
	puts "$t\t[lindex $p 2]\t[lindex $p0 2]\t$z"
	
	set sum [expr $sum + $z]
	set sumSq [expr $sumSq + $z*$z]
	incr n
}
close $out

# Compute the mean and standard error.
set mean [expr $sum/$n]
set meanSq [expr $sumSq/$n]
set se [expr sqrt(($meanSq - $mean*$mean)/$n)]

puts ""
puts "********Results: "
puts "mean dipole: $mean"
puts "standard error: $se"
mol delete top
mol delete top



