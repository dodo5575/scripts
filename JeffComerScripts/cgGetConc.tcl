# Get the system size using measure minmax.
# to use: vmd -dispdev text -e getSize.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set conc 1.0
set selTextList {"name SOD" "name CL"}
set waterText "name H2O"
# Input:
set psf cg_ions_${conc}M.psf
set pdb cg_ions_${conc}M.pdb

mol load psf $psf pdb $pdb
set waterSel [atomselect top $waterText]
set waterN [$waterSel num]
set waterCount [expr 4*$waterN]
puts "Number of water beads: $waterN"
puts "Number of waters (hypothetically): [expr $waterCount]"
$waterSel delete

foreach s $selTextList {
    set sel [atomselect top $s]
    set n [$sel num]
    puts "Number of $s: $n"
    puts "Concentration of $s: [expr 55.523*$n/$waterCount] mol/kg"
    $sel delete
}

set all [atomselect top all]
set q [measure sumweights $all weight charge]
puts "\nTotal charge: $q"
$all delete

mol delete top
exit



