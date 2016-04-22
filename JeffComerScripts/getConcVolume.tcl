# Get the system size using measure minmax.
# to use: vmd -dispdev text -e getSize.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set selTextList {"name POT" "name CLA"}
set waterText "name OH2 POT CLA"
# Input:
set psf oligo_ions.psf
set coor oligo_ions.pdb

mol load psf $psf
mol addfile $coor
set waterSel [atomselect top $waterText]
set waterN [$waterSel num]
puts "Number of waters: $waterN"
$waterSel delete

foreach s $selTextList {
    set sel [atomselect top $s]
    set n [$sel num]
    puts "Number of $s: $n"
    puts "Concentration of $s: [expr 55.523*$n/$waterN] mol/kg"
    $sel delete
}

set all [atomselect top all]
set q [measure sumweights $all weight charge]
puts "\nTotal charge: $q"
$all delete

mol delete top
exit



