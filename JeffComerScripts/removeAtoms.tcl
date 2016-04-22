# This script will remove all residues in the selection
# from psf and pdf files.
# Use with: vmd -dispdev text -e removeResidues.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set selText "not (abs(z) < 35 and x^2+y^2 > (30 + 2*(41-30)/102*abs(z))^2)"

#Input:
set psf  silica_block.psf
set pdb  silica_block.pdb
#Output:
set finalpsf silica_block1.psf
set finalpdb silica_block1.pdb

# Obtain the {segid resid name} for the selection.
mol load psf $psf pdb $pdb
set sel [atomselect top $selText]
set atomList [lsort -unique [$sel get {segname resid name}]]
set nAtoms [$sel num]
set nResidues [llength $atomList]
$sel delete

package require psfgen 1.3
resetpsf

readpsf $psf
coordpdb $pdb

# Delete the selection.
foreach atom $atomList {
	delatom [lindex $atom 0] [lindex $atom 1] [lindex $atom 2]
}

writepsf $finalpsf
writepdb $finalpdb
puts ""
puts "$nAtoms atoms were deleted."
puts "$nResidues residues were deleted."
mol delete top
exit



