# Use with: vmd -dispdev text -e generateHairpin.tcl
# Author: Jeff Comer <jcomer2@uiuc.edu>

set maxTries 100
set seg AMI
set baseName SI1
set endName N1
set range 2
# Input:
set psf amine_coating1.psf
set pdb amine_coating1.pdb
# Output:
set outName amine_coating_clean

source vector.tcl

mol load psf $psf pdb $pdb
set sel [atomselect top "segname $seg and name $baseName"]
foreach zero {0} {
    set resList [$sel get resid]
    set posList [$sel get {x y z}]
}
$sel delete

foreach res $resList pos $posList {
    set try 0
    set sel [atomselect top "segname $seg and resid $res and name $endName and within $range of not (segname $seg and resid $res)"]
    
    while {[$sel num] > 0} {
	set s [atomselect top "segname $seg and resid $res"]
	set rot [matMake4 [matRandomRot]]
	$s moveby [vecinvert $pos]
	$s move $rot
	$s moveby $pos
	$s delete
 
	$sel delete
	set sel [atomselect top "segname $seg and resid $res and name $endName and within $range of not (segname $seg and resid $res)"]
	incr try
	if {$try >= $maxTries} {
	    puts "Could not resolve clash, residue $res"
	    break
	}
    }
}

set all [atomselect top all]
$all writepdb $outName.pdb
$all writepsf $outName.psf

exit
