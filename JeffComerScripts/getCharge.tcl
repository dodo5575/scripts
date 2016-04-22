# Get the system size using measure minmax.
# to use: vmd -dispdev text -e getSize.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set selText "all"
# Input:
set psf surf_dna.psf
set pdb surf_dna.pdb


mol load psf $psf pdb $pdb
set sel [atomselect top $selText]
set q [measure sumweights $sel weight charge]
puts "\nTotal charge: $q"

$sel delete
mol delete top
exit



