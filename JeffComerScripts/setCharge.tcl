# Use with: vmd -dispdev text -e setCharge.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set selList {"type SI" "type OSI"}
set chargeList {0.9 -0.45}

#Input:
set psf silica_pore.psf
set pdb silica_pore.pdb
#Output:
set outName silica_pore_charge

# Obtain the {segid resid name} for the selection.
mol load psf $psf pdb $pdb

foreach sel $selList charge $chargeList {
    set s [atomselect top $sel]
    $s set charge $charge
    puts "Set the charge of [$s num] atoms defined by $sel to $charge."
    $s delete
}

set all [atomselect top all]
$all writepsf $outName.psf
$all writepdb $outName.pdb

mol delete top
exit
