# This script will remove water from psf and pdf files.
# Use with: vmd -dispdev text -e removeWater.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set selText "abs(x) > 20 or abs(y) > 20 or abs(z-55) > 20"

#Input:
set psf  1p2y.psf
set pdb  1p2y.pdb
set coords eq4.restart.coor
#Output:
set finalpsf soln.psf
set finalpdb soln.pdb

# Obtain the {segid resid name} for the selection.
mol load psf $psf pdb $pdb
mol addfile $coords waitfor all
set sel [atomselect top $selText]
set atomList [lsort -unique [$sel get {segname resid}]]

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



