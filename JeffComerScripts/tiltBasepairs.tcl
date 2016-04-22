f# Author: Jeff Comer <jcomer2@illinois.edu>

set structList {pore_at_basepair pore_gc_basepair}
set outList {five_prime_at five_prime_gc}
set angle 50.0; # Angle base normal makes with z axis
set selText "segname ADNA BDNA"

foreach structName $structList outName $outList {

    # Put the basepairs in the correct positions.
    mol load psf ${structName}.psf pdb ${structName}.pdb
    set sel [atomselect top $selText]
    $sel move [transaxis x [expr {-2.0*$angle}]]
    $sel delete

    set all [atomselect top all]
    $all writepsf $outName.psf
    $all writepdb $outName.pdb
    mol delete top
}
exit

