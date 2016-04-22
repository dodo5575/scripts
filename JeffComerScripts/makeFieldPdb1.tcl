# This script sets the beta value to 1.0 if force should be exerted, otherwise 0.0.
# It also puts the charge in the occupancy column.
# use with: vmd -dispdev text -e makeFieldPdb.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set seltext "protein or ions"

#Input:
set pdb  1p2y_reduced.pdb
set psf  1p2y_reduced.psf
#Output:
set finalpdb fieldIons.pdb

mol load psf $psf pdb $pdb

# Set the beta of all to 0.0 and the occupancy to the charge.
set all [atomselect top "all"]
$all set beta 0.0
$all set occupancy [$all get charge]

puts \
"The occupancy column of the pdb contains the charge gleaned from the psf."

# Set the beta of the protein to 1.0.
set sel [atomselect top $seltext]
$sel set beta 1.0

$all writepdb $finalpdb
exit



