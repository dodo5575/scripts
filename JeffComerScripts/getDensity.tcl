# Get the system size using measure minmax.
# to use: vmd -dispdev text -e getSize.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set x0 20
set x1 35
set y0 -35
set y1 35
set z0 -20
set z1 20
set selText "resname SIO2"
# Input:
set psf trap_periodic.psf
set coord trap_periodic.pdb

set dalton 1.66054e-27; # kg

mol load psf $psf
mol addfile $coord
set sel [atomselect top "($selText) and x > $x0 and x < $x1 and y > $y0 and y < $y1 and z > $z0 and z < $z1"]
puts "Number: [$sel num]"
set volume [expr ($x1-$x0)*($y1-$y0)*($z1-$z0)]
puts "Volume: $volume A^3"
set mass [measure sumweights $sel weight mass]
puts "Mass: $mass Da"
puts "Density: [expr $mass/$volume] Da/A^3"
puts "Density: [expr 1e3*$mass*$dalton/($volume*1e-24)] g/cm^3"

$sel delete
mol delete top
exit



