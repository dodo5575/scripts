# This script will remove all residues in the selection
# from psf and pdf files.
# Use with: vmd -dispdev text -e removeResidues.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set name anneal
set selText "segname V0"
set selText0 "segname U0"
set dz -0.05
set rad 1.0

#Input:
set psf combine_${name}.psf 
set pdb combine_${name}.pdb
#Output:
set finalpsf flip_${name}.psf 
set finalpdb flip_${name}.pdb

# Obtain the {segid resid name} for the selection.
mol load psf $psf pdb $pdb

set sel [atomselect top $selText]
set cen [measure center $sel weight mass]
$sel moveby [vecinvert $cen]
$sel move [transaxis y 180]

set sel0 [atomselect top $selText0]

set violate [atomselect top "($selText) and within $rad of ($selText0)"]
set num [$violate num]
$violate delete
while {$num > 0} {
    $sel moveby [list 0 0 $dz]
    
    set violate [atomselect top "($selText) and within $rad of ($selText0)"]
    set num [$violate num]
    $violate delete
}


set all [atomselect top all]
set cen [measure center $all weight mass]
$all moveby [vecinvert $cen]
$all writepsf $finalpsf
$all writepdb $finalpdb
$all delete

$sel0 delete
$sel delete
mol delete top
exit
