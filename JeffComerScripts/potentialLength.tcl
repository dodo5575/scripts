# externalForce.tcl
# Add an external electric field to a potential grid.
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
set lz [lindex $extent 2]
set dz [lindex $voxel 2]

set profileFile ${prefix}_z.dat
puts "Generating profile along z-axis: $profileFile"
set z0 [expr [lindex $g(origin) 2] + 2.0*$dz]
set z1 [expr [lindex $g(origin) 2] + $lz - 2.0*$dz]
potentialProfile g [list 0 0 $z0] [list 0 0 $z1] 100 $profileFile

puts "Length along z-axis: $lz"
exit



