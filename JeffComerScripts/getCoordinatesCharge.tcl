# Extract the coordinates of the chosen atoms.
# Use with: vmd -dispdev text -e getCoordinates.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

foreach nuc {ade thy gua cyt} {
    set selText "segname ADNA BDNA"
    #Input:
    set psf nucleo_${nuc}_pot.psf
    set coords nucleo_${nuc}_pot.pdb
    #Output:
    set outFile nucleo_${nuc}_charge.dat

    mol load psf $psf
    mol addfile $coords
    set sel [atomselect top $selText]
    foreach zero {0} {set pos [$sel get {charge x y z}]}
    set n [$sel num]
    $sel delete
    mol delete top

    set out [open $outFile w]
    foreach r $pos {
	foreach {q x y z} $r {break}
	puts $out "$q $x $y $z"
    }
    close $out

    puts "Wrote $n coordinates to $outFile."
}
exit
