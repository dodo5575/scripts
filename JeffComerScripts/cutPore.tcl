# This script will remove all residues in the selection
# from psf and pdf files.
# Use with: vmd -dispdev text -e removeResidues.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set angle 10.0
set constrict 24.0; # in angstroms

#Input:
set psf block_d8_h27_types.psf
set pdb block_d8_h27.pdb
#Output:
set finalpsf pore2.0_long.psf
set finalpdb pore2.0_long.pdb

set pi [expr 4.0*atan(1.0)]
set gamma [expr $angle*$pi/180.0]
set m [expr tan($gamma)]
set selText "sqrt(x^2 + y^2) < 0.5*$constrict + abs(z)*$m"

# Obtain the {segid resid name} for the selection.
mol load psf $psf pdb $pdb
set sel [atomselect top $selText]
foreach null {0} {set atomList [lsort -unique [$sel get {segname resid name}]]}
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



