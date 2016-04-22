# externalForce.tcl
# Add an external electric field to a potential grid.
# Requires gridForce.tcl.
# Author: Jeff Comer <jcomer2@illinois.edu>

if {$argc != 2} {
	puts "This script requires two arguments: "
	puts "  file name prefix (i.e., without .dx)"
	puts "  voltage drop"
	puts ""
	exit
}

proc trimExtension {name} {
    set ind [string last "." $name]
    return [string range $name 0 [expr $ind-1]]
}

# Parameters:
set prefix [trimExtension [lindex $argv 0]]
set voltage [lindex $argv 1]

source $env(HOME)/scripts/vector.tcl
source $env(HOME)/scripts/gridForce.tcl

set inFile ${prefix}.dx
puts "Reading $inFile"

readDx g $inFile
puts "origin: $g(origin)"
puts "number of nodes: [llength $g(data)]"

set extent [getGridSize g]
set voxel [getVoxelSize g]
set lz [lindex $extent 2]
set dz [lindex $voxel 2]

puts "Shifting the origin to zero potential."
#set v0 [averageCrossSection g 0]
set v0 [lindex $g(data) end]
addConstantPotential g [expr -$v0]

set profileFile ${prefix}_z.dat
puts "Generating profile along z-axis: $profileFile."
set z0 [expr [lindex $g(origin) 2] + 2.0*$dz]
set z1 [expr [lindex $g(origin) 2] + $lz - 2.0*$dz]
potentialProfile g [list 0 0 $z0] [list 0 0 $z1] 100 $profileFile

puts "Adding a potential drop of $voltage along the z-axis (length $lz)."
set eField0 [list 0 0 [expr $voltage/$lz]]
addConstantForce g $eField0

set profileFile1 ${prefix}_ext_z.dat
puts "Generating a new profile along z-axis: $profileFile1."
potentialProfile g [list 0 0 $z0] [list 0 0 $z1] 100 $profileFile1

writeDx g ${prefix}_ext.dx
exit



