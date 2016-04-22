# Author: Jeff Comer <jcomer2@illinois.edu>

# Parameters:
set dl 2
set ax {60 0 0}
set ay {0 60 0}
set az {0 0 138.4}
# Output:
set outName anthony.dx

source $env(HOME)/scripts/vector.tcl
source $env(HOME)/scripts/gridForce.tcl

newGridBox grid $ax $ay $az $dl
writeDx grid $outName
