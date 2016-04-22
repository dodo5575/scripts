# Extract the coordinates of the chosen atoms.
# Use with: vmd -dispdev text -e getCoordinates.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set z0 44
set rad 18
set selText "resname SIN and abs(z) > $z0 and x^2+y^2>${rad}^2"

#Input:
set psf pore1.2_all.psf
set coords output/coil_p1.2_4Va9.restart.coor
#Output:
set outFile pore1.2_surf_coords.txt

mol load psf $psf
mol addfile $coords
set sel [atomselect top $selText]
foreach zero {0} {set pos [$sel get {x y z}]}
set n [$sel num]
$sel delete
mol delete top

set out [open $outFile w]
foreach r $pos {
    foreach {x y z} $r {break}
    puts $out "$x $y $z"
}
close $out

puts "Wrote $n coordinates to $outFile."
exit



