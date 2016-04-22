# Author: Jeff Comer <jcomer2@illinois.edu>

# Parameters:
set dl 1.0
set lx 12
set ly 12
set lz 12
# Output:
set outName diffGrid.dx

source $env(HOME)/scripts/vector.tcl
source $env(HOME)/scripts/gridForce.tcl

newGridFit grid $lx $ly $lz $dl
writeDx grid $outName
