# Compute the solvent-accessible surface area for a trajectory.
# to use: vmd -dispdev text -e sasaDcd.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

# SASA radius
set radius 2.0
set dcdFreq 5000
set selText "protein"
set startFrame 0
set timestep 1.0
set stride 10

# Input:
set pdb  1p2y_reduced.pdb
set psf  1p2y_reduced.psf
set dcdPrefix null
set dcdSuffix ".dcd"
set dcdSet {0 1 2 3}
# Output:
set outFile sasa_no_null.txt

# Get the time change between frames in femtoseconds.
set dt [expr $timestep*$dcdFreq]
set dcdEnd [llength $dcdSet]

# Open the output file.
set out [open $outFile w]

# Load the system.
mol load psf $psf pdb $pdb
set sel [atomselect top $selText]

# Loop over the dcd files.
set nFrames0 0
for {set dcd 0} {$dcd < $dcdEnd} {incr dcd} {

# Load the trajectory.
animate delete all
set dcdNum [lindex $dcdSet $dcd]
mol addfile "${dcdPrefix}${dcdNum}${dcdSuffix}" waitfor all
set nFrames [molinfo top get numframes]
puts [format "Reading %i frames." $nFrames]

# Start at "startFrame" and move forward, computing
# the sasa at each step.
for {set f $startFrame} {$f < $nFrames} {incr f $stride} {
	molinfo top set frame $f

	set a [measure sasa $radius $sel]
		
	# Get the time in nanoseconds for this frame.
	set t [expr ($nFrames0+$f+1)*$dt*1.e-6]
	
	puts $out "$t $a"
	puts -nonewline [format "FRAME %i: " $f]
	puts "$t $a"
	
}
set nFrames0 [expr $nFrames+$nFrames0]
}

close $out
exit




