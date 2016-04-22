# This script will remove water from psf and pdf files.
# Use with: vmd -e drawPore.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

# Parameters:
set radius0 8.0
set radius1 16.8
set length 100.0
set mat Opaque
set col cyan

source drawPore.tcl
set molid top

graphics $molid delete all
graphics $molid material $mat
graphics $molid color $col
drawPore top $radius0 $radius1 $length 40



