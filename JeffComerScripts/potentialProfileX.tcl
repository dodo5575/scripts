# Requires gridForce.tcl.
# Author: Jeff Comer <jcomer2@illinois.edu>

if {$argc != 1} {
	puts "This script requires two arguments: "
	puts "  file name prefix (i.e., without .dx)"
	puts ""
	exit
}

proc trimExtension {name} {
    set ind [string last "." $name]
    return [string range $name 0 [expr $ind-1]]
}
set prefix [trimExtension [lindex $argv 0]]

source $env(HOME)/scripts/vector.tcl
source $env(HOME)/scripts/gridForce.tcl

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

set profileFile ${prefix}_x.dat
puts "Generating profile along x-axis: $profileFile"

set y [expr [lindex $g(origin) 1] + 0.5*$ly]
set z [expr [lindex $g(origin) 2] + 0.5*$lz]
puts "Position: x $y $z"

set x0 [expr [lindex $g(origin) 0] + 2.0*$dx]
set x1 [expr [lindex $g(origin) 0] + $lx - 2.0*$dx]
potentialProfile g [list $x0 $y $z] [list $x1 $y $z] 100 $profileFile
exit


