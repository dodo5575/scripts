# This script will remove water from psf and pdf files.
# Use with: vmd -e drawPore.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

# Parameters:
set diameter0 24.15
set diameter1 61.78
set delta 1.0
set length 200.0
set mat Translucent
set col cyan

source drawPore.tcl

set l [expr $length+$delta]
set d0 [expr 0.5*$diameter0-$delta]
set d1 [expr 0.5*$diameter1-$delta]
set molid top

graphics $molid delete all
graphics $molid material $mat
graphics $molid color $col
drawPore top $d0 $d1 $l 40



