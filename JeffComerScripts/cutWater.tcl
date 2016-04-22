# This script will remove water from psf and pdf files.
# Use with: vmd -dispdev text -e cutWater.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set selText "abs(z) > 55 and (water or ions)"

#Input:
set psf  1p2y.psf
set pdb  1p2y.pdb
set bincoords eq4.restart.coor
#Output:
set finalpsf 1p2y_reduced.psf
set finalpdb 1p2y_reduced.pdb

# Obtain the {segid resid} for the selection.
mol load psf $psf pdb $pdb
mol addfile $bincoords waitfor all
set sel [atomselect top $selText]
set atomList [lsort -unique [$sel get {segid resid}]]

package require psfgen 1.3
resetpsf

readpsf $psf
coordpdb $pdb

# Delete the selection.
foreach atom $atomList {
	delatom [lindex $atom 0] [lindex $atom 1]
}

writepsf $finalpsf
writepdb $finalpdb
exit



