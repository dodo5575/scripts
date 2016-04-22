# potentialVmdToVolts.tcl
# Convert from k_B*T/e (300 K) to volts.
# Requires gridForce.tcl.
# Author: Jeff Comer <jcomer2@illinois.edu>

if {$argc != 1} {
	puts "This script requires one argument: "
	puts "  file name"
	puts ""
	exit
}

proc trimExtension {name} {
    set ind [string last "." $name]
    return [string range $name 0 [expr $ind-1]]
}

# Parameters:
set prefix [trimExtension [lindex $argv 0]]

source $env(HOME)/scripts/vector.tcl
source $env(HOME)/scripts/gridForce.tcl

set inFile ${prefix}.dx
puts "Reading $inFile"

readDx g $inFile
scalePotential g 0.0258520296
writeDx g ${prefix}_volts.dx
