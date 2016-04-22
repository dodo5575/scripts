# Add harmonic constraints to silicon nitride.
# to use: vmd -dispdev text -e langevinWater.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

# Parameters:
# Langevin constant in ps^-1.
set dampingList {108.05202 101.77514}
set selTextList {"name POT" "name CLA"}
set prefix 100mM_100Da
set outName 100mM_langevin

set psf $prefix.psf
set pdb $prefix.pdb
mol load psf $psf pdb $pdb
set selAll [atomselect top all]

# Set the spring constants to zero for all atoms.
$selAll set occupancy 0.0
$selAll set beta 0.0

foreach selText $selTextList damping $dampingList {
    set sel [atomselect top $selText]
    $sel set beta $damping

    puts "Set beta to $damping for [$sel num] atoms given by \"$selText\"."
    $sel delete
}

$selAll writepdb $outName.pdb
$selAll delete
mol delete top
exit



