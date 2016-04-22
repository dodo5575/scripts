# potentialGhost.tcl
# Add ghost points.
# Requires gridForce.tcl.
# Author: Jeff Comer <jcomer2@illinois.edu>

if {$argc != 1} {
	puts "This script requires one arguments: "
	puts "  file name prefix (i.e., without .dx)"
        puts ""
	exit
}

# Parameters:
set prefix [lindex $argv 0]

source gridForce.tcl

set inFile ${prefix}.dx
puts "Reading $inFile"

readDx g $inFile
ghostGrid g gg
writeDx gg ${prefix}_ghost.dx



