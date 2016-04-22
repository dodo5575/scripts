# This script will remove all residues in the selection
# from psf and pdf files.
# Use with: vmd -dispdev text -e removeResidues.tcl
# Author: Jeff Comer <jeffcomer at gmail>

set selText "not segname DEN"

#Input:
set psf dna1537+pamam_eq2.psf
set pdb dna1537+pamam_eq2.pdb
#Output:
set outName pamam

set finalPsf $outName.psf 
set finalPdb $outName.pdb

# Obtain the {segid resid name} for the selection.
mol load psf $psf pdb $pdb
set sel [atomselect top $selText]
foreach zero {0} {set atomList [lsort -unique [$sel get {segname resid}]]}
set nAtoms [$sel num]
set nResidues [llength $atomList]
$sel delete

package require psfgen 1.3
resetpsf

readpsf $psf
coordpdb $pdb

# Delete the selection.
foreach atom $atomList {
	delatom [lindex $atom 0] [lindex $atom 1]
}

writepsf $finalPsf
writepdb $finalPdb
puts ""
puts "$nAtoms atoms were deleted."
puts "$nResidues residues were deleted."
mol delete top
exit




