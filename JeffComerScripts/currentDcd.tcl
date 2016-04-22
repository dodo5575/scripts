# Calculate the current for a trajectory.
# Writes time(ns) current(nA) to a file.
# to use: vmd -dispdev text -e currentDcd.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set dcdFreq 5000
set selText "ions"
set timestep 1.0
set startTime 1.277

# Input:
set pdb pore+dna-all.pdb
set psf pore+dna-all.psf
set dcdPrefix run
set dcdSuffix "_4V.dcd"
set dcdSet {3 4 5 6 7}
set xsc run6_4V.restart.xsc
# Output:
set outFile curr_conform2_4V_[lindex $dcdSet 0]-[lindex $dcdSet end].txt

set dcdEnd [llength $dcdSet]
# Get the time change between frames in femtoseconds.
set dt [expr $timestep*$dcdFreq]

# Read the system size from the xsc file.
# Note: This only works for lattice vectors along the axes!
set in [open $xsc r]
foreach line [split [read $in] "\n"] {
	if {![string match "#*" $line]} {
		set param [split $line]
		puts $param
		set lx [lindex $param 1]
		set ly [lindex $param 5]
		set lz [lindex $param 9]
		break
	}
}
puts "NOTE: The system size is $lx $ly $lz.\n"
close $in

# Open the output file.
set out [open $outFile w]
#puts $out "sum of q*v for $psf with trajectory $dcd"
#puts $out "t(ns) I(A)"

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

for {set i 0} {$i < 1} {incr i} {
	# Get the charge of each atom.
	set q [$sel get charge]

	# Get the position data for the first frame.
	molinfo top set frame 0
	set z0 [$sel get z]
}

# Start at "startFrame" and move forward, computing
# current at each step.
for {set f 1} {$f < $nFrames} {incr f} {
	molinfo top set frame $f
	
	# Get the position data for the current frame.
	set z1 [$sel get z]
	
	# Find the displacements in the z-direction.
	set dz {}
	foreach r0 $z0 r1 $z1 {
		# Compensate for jumps across the periodic cell.
		set z [expr $r1-$r0]
		if {[expr $z > 0.5*$lz]} {set z [expr $z-$lz]}
		if {[expr $z <-0.5*$lz]} {set z [expr $z+$lz]}
		
		lappend dz $z
	}
	
	# Compute the average charge*velocity between the two frames.
	set qvsum [expr [vecdot $dz $q] / $dt]
		
	# We first scale by the system size to obtain the z-current in e/fs.
	set currentZ [expr $qvsum/$lz]
	# Now we convert to nanoamperes.
	set currentZ [expr $currentZ*1.60217733e5]
	# Get the time in nanoseconds for this frame.
	set t [expr ($nFrames0+$f+0.5)*$dt*1.e-6 + $startTime]
			
	# Write the current.
	puts $out "$t $currentZ"
	puts -nonewline [format "FRAME %i: " $f]
	puts "$t $currentZ"
	
	# Store the position data for the next computation.
	set z0 $z1
}
set nFrames0 [expr $nFrames+$nFrames0]
}

close $out
exit




