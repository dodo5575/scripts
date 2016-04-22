# Make a pdb using NAMD restart coordinates.
# Use with: vmd -dispdev text -e makeCoordPdb.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set seltext "all"

#Input:
set psf pore_at_basepair_1M.psf
set coords output/pore_at_1M_eq5.restart.coor
#Output:
set finalpdb pore_at_1M_eq5.pdb

mol load psf $psf
mol addfile $coords waitfor all
set sel [atomselect top $seltext]
$sel writepdb $finalpdb
exit



