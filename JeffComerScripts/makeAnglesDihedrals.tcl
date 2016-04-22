# This script will generate angles and dihedrals.
# Use with: vmd -dispdev text -e makeAnglesDihedrals.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set selText "all"

#Input:
set psf  silica_pore2.psf
set pdb  silica_pore1.pdb
#Output:
set finalpsf silica_pore3.psf
set finalpdb silica_pore3.pdb

package require psfgen 1.3
resetpsf

readpsf $psf
coordpdb $pdb

regenerate angles dihedrals

writepsf $finalpsf
writepdb $finalpdb
exit



