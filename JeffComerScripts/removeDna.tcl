# This script will remove DNA from psf and pdf files.
# use with: vmd -dispdev text -e removeDNA.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set seltext "segname ADNA"

#Input:
set psf  pore+dna-all.psf
set pdb  pore+dna-all.pdb
set bincoords eq3.restart.coor
#Output:
set finalpsf pore-all.psf
set finalpdb pore-all.pdb

set temppdb tmp.pdb

# Load the molecule with the given coordinate file.
mol load psf $psf pdb $pdb
mol addfile $bincoords waitfor all
set all [atomselect top all]
$all writepdb $temppdb

# Obtain the {segid resid name} for the selection.
set sel [atomselect top $seltext]
set atomList [$sel get {segid resid name}]

package require psfgen 1.3
resetpsf

readpsf $psf
coordpdb $temppdb

# Delete the selection.
foreach atom $atomList {
	delatom [lindex $atom 0] [lindex $atom 1] [lindex $atom 2]
}

writepsf $finalpsf
writepdb $finalpdb
exit



