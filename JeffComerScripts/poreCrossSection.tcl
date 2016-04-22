# This script will remove all residues in the selection
# from psf and pdf files.
# Use with: vmd -dispdev text -e removeResidues.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

#Input:
set pore 5
set buffFrac 0.3
set probe 0.0; # rMin/2 = 1.7683 for TIP3P oxygen
set nameList {pore1.2 pore1.6 pore2.0 pore2.0_thick pore2.0_trap pore2.2}
set name [lindex $nameList $pore]
set pi [expr 4.0*atan(1.0)]

set inFile dist_${name}.dx
set outArea [open area_${name}.dat w]
set outRad [open diam_${name}.dat w]

source vector.tcl
source gridForce.tcl

# Initialize the grid.
readDx grid $inFile
puts "Read $inFile."
set nx $grid(nx)
set ny $grid(ny)
set nz $grid(nz)
puts "Grid dimensions are $grid(nx) $grid(ny) $grid(nz)."

# Select the region to be scanned.
set ix0 [expr int($buffFrac*$nx)]
set ix1 [expr $nx - $ix0]
set iy0 [expr int($buffFrac*$ny)]
set iy1 [expr $ny - $iy0]
set iz0 [expr int($buffFrac*$nz)]
set iz1 [expr $nz - $iz0]

set z0 [expr 0.1*[lindex $grid(origin) 2]]
puts "Origin along z: $z0 nm"
set dz [expr 0.1*[lindex $grid(delta) 2 2]]
puts "Grid spacing along z: $dz nm"
set normal [vecCross [lindex $grid(delta) 0] [lindex $grid(delta) 1]]
set voxelArea [expr 0.01*abs([vecLength $normal])]
puts "Voxel xy-area: $voxelArea nm^2"

for {set iz $iz0} {$iz < $iz1} {incr iz} {
    set count 0
    
    for {set iy $iy0} {$iy < $iy1} {incr iy} {
	for {set ix $ix0} {$ix < $ix1} {incr ix} {
	    set i [expr $iz + $grid(nz)*$iy + $grid(nz)*$grid(ny)*$ix]
	    set d [lindex $grid(data) $i]

	    if {$d > $probe} {
		incr count
	    }
	}
    }

    set z [expr $z0 + $iz*$dz]
    set area [expr $count*$voxelArea]
    set diam [expr 2.0*sqrt($area/$pi)]

    puts "$z $area $diam"
    puts $outArea "$z $area"
    puts $outRad "$z $diam"
}
close $outArea
close $outRad
exit



