# This script will remove all residues in the selection
# from psf and pdf files.
# Use with: vmd -dispdev text -e cutRectangle.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set selText "resname SIO"
set solText "water or ions"

#Input:
set psf  pore-solv.psf
set pdb  pore-solv.pdb
#Output:
set finalpsf pore-rect.psf
set finalpdb pore-rect.pdb

mol load psf $psf pdb $pdb
set sel [atomselect top $selText]

set minmax [measure minmax $sel]
set x0 [lindex $minmax 0 0]
set y0 [lindex $minmax 0 1]
set x1 [lindex $minmax 1 0]
set y1 [lindex $minmax 1 1]
$sel delete

set violatorText "(${solText}) and (x < $x0 or x > $x1 or y < $y0 or y > $y1)"
set violators [atomselect top $violatorText]

set atomList [lsort -unique [$violators get {segname resid}]]
set nAtoms [$violators num]
set nResidues [llength $atomList]
$violators delete

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
puts ""
puts "$nAtoms atoms were deleted."
puts "$nResidues residues were deleted."
mol delete top
exit



