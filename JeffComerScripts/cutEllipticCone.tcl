# This script will cut an elliptical double cone from the psf and pdb.
# Use with: vmd -dispdev text -e cutEllipticCone.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

# Parameters:
set selText "segname SIN"
set slope 0.25
# x and y diameters in nanometers:
set dmin_x 3.0
set dmin_y 1.5
# Input:
set psf  pore1.5.psf
set pdb  pore1.5.pdb
# Output:
set finalpsf pore_ellipse.psf
set finalpdb pore_ellipse.pdb

set rx [expr 5.0*$dmin_x]
set ry [expr 5.0*$dmin_y]

# Obtain the {segid resid name} for the selection.
mol load psf $psf pdb $pdb
set vioText "($selText) and (x/($rx+$slope*abs(z)))^2 + \
(y/($ry+$slope*abs(z)))^2 < 1."
set sel [atomselect top $vioText]
set atomList [$sel get {segname resid name}]
set nViolators [$sel num]

puts "$nViolators atom(s) will be removed."

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
exit



