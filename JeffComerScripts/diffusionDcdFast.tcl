# Find binding events and write their durations. Do this for several radii.
# to use: vmd -dispdev text -e multicurrentDCD.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set dcdFreq 5000
set timestep 1.0
set selText "name POT"
set substrateText "resname SIN"
set window 20.

# Input:
set pdb  DS3-noDNA.pdb
set psf  DS3-noDNA.psf
set dcdPrefix  /Scr/alek/NO_DNA/output/E0.1-DS3-noDNA-
set dcdSuffix ".dcd"
set dcdSet {5}
# Output:
set outFile diffuse_DS3_fast.txt

set dcdEnd [llength $dcdSet]
# Get the time change between frames in femtoseconds.
set dt [expr $timestep*$dcdFreq]

# Find the nearest an atom in selection "sel" comes to position "r".
proc findNearestDist {r selText window} {
	# Check only those atoms within the window.
	foreach {x y z} $r {break}
	set w2 [expr $window*$window]
	set sel [atomselect top "$selText and (x-$x)^2+(y-$y)^2+(z-$z)^2 < $w2"]
	
	# Check for a quick out.
	if {[$sel num] == 0} {return [expr 2.0*$window]}
	
	set rSel [$sel get {x y z}]
	set nearDistSq [veclength2 [vecsub $r [lindex $rSel 0]]]
	foreach r0 $rSel {
		set distSq [veclength2 [vecsub $r $r0]]
		if {$distSq < $nearDistSq} {
			set nearDistSq $distSq
		}
	}
	
	$sel delete
	return [expr sqrt($nearDistSq)]
}

# Load the system.
mol load psf $psf pdb $pdb
set sel [atomselect top $selText]
set out [open $outFile w]

# Loop over the dcd files.
set nFrames0 0
for {set dcd 0} {$dcd < $dcdEnd} {incr dcd} {
puts $dcd

# Load the trajectory.
animate delete all
set dcdNum [lindex $dcdSet $dcd]
mol addfile "${dcdPrefix}${dcdNum}${dcdSuffix}" waitfor all
set nFrames [molinfo top get numframes]
puts [format "Reading %i frames." $nFrames]

foreach i {1} {
	# Get the position data for the first frame.
	molinfo top set frame 0
	set r0 [$sel get {x y z}]
}

# Move forward, computing
# current at each step.
for {set f 1} {$f < $nFrames} {incr f} {
	molinfo top set frame $f
	
	# Get the position data for the current frame.
	set r1 [$sel get {x y z}]
	
	foreach p0 $r0 p1 $r1 {
		set v [expr [veclength [vecsub $p1 $p0]]/$dt]
		set p [vecscale 0.5 [vecadd $p1 $p0]]
		set dSubstrate [findNearestDist $p $substrateText $window]
		if {$dSubstrate <= $window} {puts $out "$v $dSubstrate"}
	}
	
	puts "FRAME $f: $v $dSubstrate"
	
	# Store these coordinates for the next iteration.
	set r0 $r1	
}
set nFrames0 [expr $nFrames+$nFrames0]
}
close $out

exit




