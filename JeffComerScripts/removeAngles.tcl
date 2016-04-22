# This script will remove all residues in the selection
# from psf and pdf files.
# Use with: vmd -dispdev text -e removeResidues.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

#Input:
set psf cg_hemo_ions1.0.psf
set pdb cg_hemo_ions1.0.pdb
#Output:
set finalPsf pore+sol.psf

mol load psf $psf pdb $pdb
set all [atomselect top all]
$all writepsf $finalPsf
$all delete

mol delete top
exit



