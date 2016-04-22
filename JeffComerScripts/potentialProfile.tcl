# Requires gridForce.tcl.
# Author: Jeff Comer <jcomer2@illinois.edu>

if {$argc != 1} {
	puts "This script requires two arguments: "
	puts "  file name prefix (i.e., without .dx)"
	puts ""
	exit
}

# Parameters:
set prefix [lindex $argv 0]

source vector.tcl
source gridForce.tcl

set inFile ${prefix}.dx
puts "Reading $inFile"

readDx g $inFile
puts "origin: $g(origin)"
puts "number of nodes: [llength $g(data)]"

set extent [getGridSize g]
set voxel [getVoxelSize g]
set lx [lindex $extent 0]
set ly [lindex $extent 1]
set lz [lindex $extent 2]
set dx [lindex $voxel 0]
set dy [lindex $voxel 1]
set dz [lindex $voxel 2]

set profileFile ${prefix}_z.dat
puts "Generating profile along z-axis: $profileFile"

set x [expr [lindex $g(origin) 0] + 0.5*$lx]
set y [expr [lindex $g(origin) 1] + 0.5*$ly]
puts "Position: $x $y z"

set z0 [expr [lindex $g(origin) 2] + 2.0*$dz]
set z1 [expr [lindex $g(origin) 2] + $lz - 2.0*$dz]
potentialProfile g [list $x $y $z0] [list $x $y $z1] 100 $profileFile
exit


