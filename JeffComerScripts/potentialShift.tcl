# potentialScale.tcl
# Scale by a factor.
# Requires gridForce.tcl.
# Author: Jeff Comer <jcomer2@illinois.edu>

if {$argc != 3} {
    puts "This script requires two arguments: "
    puts "  file name"
    puts "  shift"
    puts "  output file"
    puts ""
    exit
}

# Parameters:
set prefix [lindex $argv 0]
set factor [lindex $argv 1]
set outFile [lindex $argv end]

source vector.tcl
source gridForce.tcl

set inFile $prefix
puts "Reading $inFile"

readDx g $inFile
addConstantPotential g $factor
writeDx g $outFile




