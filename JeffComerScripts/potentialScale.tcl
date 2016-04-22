# potentialScale.tcl
# Scale by a factor.
# Requires gridForce.tcl.
# Author: Jeff Comer <jcomer2@illinois.edu>

if {$argc != 2} {
	puts "This script requires two arguments: "
	puts "  file name"
	puts "  factor by which to scale"
        puts ""
	exit
}

proc trimExtension {name} {
    set ind [string last "." $name]
    return [string range $name 0 [expr $ind-1]]
}

# Parameters:
set prefix [trimExtension [lindex $argv 0]]
set factor [lindex $argv 1]

source $env(HOME)/scripts/vector.tcl
source $env(HOME)/scripts/gridForce.tcl

set inFile ${prefix}.dx
puts "Reading $inFile"

readDx g $inFile
scalePotential g $factor
writeDx g ${prefix}_scaled${factor}.dx



