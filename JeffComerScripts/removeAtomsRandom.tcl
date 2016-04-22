# Author: Jeff Comer <jcomer2@illinois.edu>

set frac 0.3
set selText "abs(z) > 10"
#Input:
set psf block_d9z8_bonds.psf
set pdb block_d9z8_bonds.pdb
#Output:
set finalpsf rough_d9z8.psf
set finalpdb rough_d9z8.pdb

proc quiet {} {
}

# Obtain the {segid resid name} for the selection.
mol load psf $psf pdb $pdb
set sel [atomselect top $selText]
set atomList [lsort -unique [$sel get {segname resid name}]]; quiet
set nAtoms [$sel num]
set nResidues [llength $atomList]
$sel delete

set violators {}
foreach atom $atomList {
    if {rand() < $frac} {
	lappend violators $atom
    }
}

package require psfgen 1.3
resetpsf

readpsf $psf
coordpdb $pdb

# Delete the selection.
foreach atom $violators {
	delatom [lindex $atom 0] [lindex $atom 1] [lindex $atom 2]
}

writepsf $finalpsf
writepdb $finalpdb
puts ""
puts "$nAtoms atoms were deleted."
puts "$nResidues residues were deleted."
mol delete top
exit
