# Calculate the current for a trajectory.
# to use: vmd -dispdev text -e poreVolume.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set nSegments 24
set binZ0 -120.0
set binZ1 120.0
set radius 2.0
set attempts 100

# Input:
set pdb  pore+dna-all.pdb
set psf  pore+dna-all.psf
set coor run2_4V.restart.coor
# Output:
set outFile vol_null100.txt

# Load the system.
mol load psf $psf pdb $pdb
mol addfile $coor waitfor all
set all [atomselect top all]
set minmax [measure minmax $all]

# Prepare the segments.
set binDz [expr ($binZ1-$binZ0)/$nSegments]
set binZ {}
for {set i 0} {$i < $nSegments} {incr i} {
	lappend binZ [expr $binZ0+$i*$binDz] [expr $binZ0+($i+1)*$binDz]
}

set out [open $outFile w]
set rad2 [expr $radius*$radius]
set x0 [lindex $minmax 0 0]
set x1 [lindex $minmax 1 0]
set y0 [lindex $minmax 0 1]
set y1 [lindex $minmax 1 1]
set dx [expr $x1-$x0]
set dy [expr $y1-$y0]
set vol0 [expr $binDz*$dx*$dy]
puts "\nComputing the volume of each segment..."

proc probe {x y z rad2} {
	set sel [atomselect top "(water or ions) and (x-$x)^2 + (y-$y)^2 + (z-$z)^2 < $rad2"]
	return [$sel num]
}

# Probe to determine the volume.
foreach {z0 z1} $binZ {
	set in 0
	for {set i 0} {$i < $attempts} {incr i} {
		set x [expr rand()*$dx + $x0]
		set y [expr rand()*$dy + $y0]
		set z [expr rand()*$binDz + $z0]

		if {[expr [probe $x $y $z $rad2] > 0]} {
			incr in
		}
	}
	
	set vol [expr $vol0*$in/floor($attempts)]
	puts $vol
	puts $out $vol
}
close $out

exit




