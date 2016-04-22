# Get the system size using measure minmax.
# to use: vmd -dispdev text -e getSize.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set desiredConc 1.0; # mol
set selTextList {"name POT" "name CLA"}
set waterText "name OH2"
# Input:
set psf pore+dna2.psf
set coor pore+dna2.pdb

mol load psf $psf
mol addfile $coor
set waterSel [atomselect top $waterText]
set waterN [$waterSel num]
puts "Number of waters: $waterN"
$waterSel delete

foreach s $selTextList {
    set sel [atomselect top $s]
    set n [$sel num]
    puts "\nNumber of $s: $n"
    puts "Desired number of $s: [expr $desiredConc*$waterN/55.523]"
    puts "Concentration of $s: [expr 55.523*$n/$waterN] mol/kg"
    $sel delete
}

set all [atomselect top all]
set q [measure sumweights $all weight charge]
puts "\nTotal charge: $q"
$all delete

mol delete top
exit



