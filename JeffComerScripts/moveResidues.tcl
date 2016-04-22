# This script combines two pdb/psf pairs.
# Use with: vmd -dispdev text -e combine.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set lx 25
set ly 25
set lz 87.0097
set selText "resname DMMP"
set waterText "water"
# Input:
set psf conc_middling.psf
set pdb conc_middling.pdb
# Output:
set outName conc_middling_fold

proc randomPos {lx ly lz} {
    set x [expr (rand()-0.5)*$lx]
    set y [expr (rand()-0.5)*$ly]
    set z [expr (rand()-0.5)*$lz]

    return [list $x $y $z]
}

# Load the structures.
mol load psf $psf pdb $pdb

# Set the periodic cell.
molinfo top set alpha 90
molinfo top set beta 90
molinfo top set gamma 90
molinfo top set a $lx
molinfo top set b $ly
molinfo top set c $lz

# Get the residues.
set sel [atomselect top $selText]
set resList [lsort -unique -integer [$sel get resid]]
$sel delete

foreach res $resList {
    set thisText "($selText) and resid $res"
    set s [atomselect top $thisText]
    
    set cen [measure center $s weight mass]
    $s moveby [vecsub [randomPos $lx $ly $lz] $cen]
    set bad [atomselect top "not (($thisText) or ($waterText)) and pbwithin 3.0 of ($thisText)"]
    set badNum [$bad num]
    $bad delete

    while {$badNum > 0} {
	set cen [measure center $s weight mass]
	$s moveby [vecsub [randomPos $lx $ly $lz] $cen]
	set bad [atomselect top "not (($thisText) or ($waterText)) and pbwithin 3.0 of ($thisText)"]
	set badNum [$bad num]
	$bad delete
    }    
}

set all [atomselect top all]
$all writepdb $outName.pdb
$all delete
mol delete top
exit
