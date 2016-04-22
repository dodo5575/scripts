# Neutralize a psf.
# Use with: vmd -dispdev text -e .tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set q 0.0
set selText all
#Input:
set psf silica_block_bonds_charged.psf
set pdb silica_block_bonds.pdb
#Output:
set finalPsf silica_block_bonds.psf

# Obtain the {segid resid name} for the selection.
mol load psf $psf pdb $pdb
set all [atomselect top all]
set sel [atomselect top $selText]
$sel set charge $q
puts "Neutralizing [$sel num] atoms defined by $selText"
puts "Warning! All angles, dihedrals, and impropers will be lost."
$all writepsf $finalPsf

$sel delete
$all delete
mol delete top
exit



