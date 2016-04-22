# Extract the coordinates of the chosen atoms.
# Use with: vmd -dispdev text -e getCoordinates.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set selText "nucleic"
#Input:
set psf run_at.psf
set coords run_at.pdb
#Output:
set outFile dna_at.txt

mol load psf $psf
mol addfile $coords
set sel [atomselect top $selText]
foreach zero {0} {set pos [$sel get {x y z radius}]}
set n [$sel num]
$sel delete
mol delete top

set out [open $outFile w]
foreach r $pos {
    foreach {x y z rad} $r {break}
    puts $out "$x $y $z $rad"
}
close $out

puts "Wrote $n coordinates to $outFile."
exit



