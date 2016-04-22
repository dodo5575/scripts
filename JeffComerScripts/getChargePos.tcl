# Extract the coordinates of the chosen atoms.
# Use with: vmd -dispdev text -e getCoordinates.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set selText "segname ADNA BDNA"

#Input:
set psf nw_dna_water67700_mem10.psf
set coords nw_dna_mem10_last.dcd
#Output:
set outName nw_dna_mem10_last

set outCharge $outName.charge.dat
set outPos $outName.pos.dat

mol load psf $psf
mol addfile $coords
set sel [atomselect top $selText]
foreach zero {0} {
    set pos [$sel get {x y z}]
    set charge [$sel get charge]
}
set n [$sel num]
$sel delete
mol delete top

set outQ [open $outCharge w]
set outR [open $outPos w]
foreach r $pos q $charge {
    foreach {x y z} $r {break}
    puts $outR "$x $y $z"
    puts $outQ "$q"
}
close $outQ
close $outR

puts "Wrote $n coordinates to $outPos and $n charges to $outCharge."
exit


