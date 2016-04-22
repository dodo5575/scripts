# Calculate the sum of charge*velocity of each atom for a trajectory.
# to use: vmd -dispdev text -e chargeVelocity.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set dcdFreq 100
set selText "ions"
set startFrame 0

# Input:
set pdb  lipid-all.pdb
set psf  lipid-all.psf
set dcd  run6.dcd
set xsc run6.restart.xsc
# Output:
set outFile run6_freq0.1.txt

# Get the time change between frames in picoseconds.
# This assumes a timestep of 1.0 fs.
set dt [expr $dcdFreq/1000.]

# Read the system size from the xsc file.
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
#puts $out "sum of q*v for $psf with trajectory $dcd"
#puts $out "t(ps) q*v(e A/ps)"

for {set i 0} {$i < 1} {incr i} {
	# Get the charge of each atom.
	set qList [$sel get charge]

	# Get the position data for the first frame.
	molinfo top set frame $startFrame
	set r0List [$sel get {x y z}]
}

# Start at "startFrame" and move forward "stride" frames, computing
# q*v at each step.
set n 1
for {set f [expr $startFrame+1]} {$f < $nFrames && $n > 0} {incr f} {
	molinfo top set frame $f
	
	# Get the position data for the current frame.
	set r1List [$sel get {x y z}]
	
	# Compute the sum of q*v over the atoms.
	set qvsum [veczero]
	foreach r0 $r0List r1 $r1List q $qList {
		# Compensate for jumps across the periodic cell.
		set dx [expr [lindex $r1 0] - [lindex $r0 0]]
		set dy [expr [lindex $r1 1] - [lindex $r0 1]]
		set dz [expr [lindex $r1 2] - [lindex $r0 2]]
		
		if {[expr $dx > 0.5*$lx]} {set dx [expr $dx-$lx]}
		if {[expr $dx <-0.5*$lx]} {set dx [expr $dx+$lx]}
		if {[expr $dy > 0.5*$ly]} {set dy [expr $dy-$ly]}
		if {[expr $dy <-0.5*$ly]} {set dy [expr $dy+$ly]}
		if {[expr $dz > 0.5*$lz]} {set dz [expr $dz-$lz]}
		if {[expr $dz <-0.5*$lz]} {set dz [expr $dz+$lz]}
		
		# Compute the average velocity between the two frames.
		set qv [vecscale [list $dx $dy $dz] [expr $q/$dt]]		
		set qvsum [vecadd $qvsum $qv]
	}
	set t [expr ($f+0.5)*$dt]
	
	# Write the q*v_z.
	puts $out "$t [lindex $qvsum 2]"
	#puts -nonewline [format "FRAME %i: " [molinfo top get frame]]
	#puts "$t $qvsum"
	
	# Store the postion data for the next computation.
	set r0List $r1List
}
close $out
exit




