# potentialVmdToVolts.tcl
# Covert from k_B*T/e (300 K) to volts.
# Requires gridForce.tcl.
# Author: Jeff Comer <jcomer2@illinois.edu>

if {$argc != 1} {
	puts "This script requires one argument: "
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
puts "number of nodes: [llength $g(data)]"

set extent [getGridSize g]
set voxel [getVoxelSize g]
set lz [lindex $extent 2]
set dz [lindex $voxel 2]

puts "system length along z-axis: $lz"

set profileFile ${prefix}_z.dat
puts "Generating profile along z-axis: $profileFile"
set z0 [expr [lindex $g(origin) 2] + 2.0*$dz]
set z1 [expr [lindex $g(origin) 2] + $lz - 2.0*$dz]
potentialProfile g [list 0 0 $z0] [list 0 0 $z1] 100 $profileFile



