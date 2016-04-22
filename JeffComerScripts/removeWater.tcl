# This script will remove water from psf and pdf files.
# Use with: vmd -dispdev text -e removeWater.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set seltext "water or ions"

#Input:
set psf  1p2y.psf
set pdb  1p2y.pdb
#Output:
set finalpsf 1p2y_vacuum.psf
set finalpdb 1p2y_vacuum.pdb

# Obtain the {segid resid name} for the selection.
mol load psf $psf pdb $pdb
set sel [atomselect top $seltext]
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



