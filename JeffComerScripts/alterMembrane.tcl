# Add harmonic constraints to silicon nitride.
# to use: vmd -dispdev text -e constrainSilicon.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

# Parameters:
set keyList {14 14.5 15 15.5 16 16.5 17 17.5 18 18.5 19}
set selText "resname SIO2"
# Input:
set memPsf DMMP_sol_raw.psf
set memPdb DMMP_sol_raw.pdb
set psf DMMP_sol.psf
set pdbPrefix wham_dmmp_z
# Output:
set outPrefix wham_raw_z

# Load the membrane coordinates.
mol load psf $memPsf pdb $memPdb
set sel [atomselect top $selText]
foreach zero {0} {set memPos [$sel get {x y z}]}
$sel delete
mol delete top


mol load psf $psf
set all [atomselect top all]
set sel [atomselect top $selText]
foreach key $keyList {
    animate delete all
    set pdb ${pdbPrefix}${key}.pdb
    mol addfile $pdb

    $sel set {x y z} $memPos
    $all writepdb ${outPrefix}${key}.pdb

}

$all delete
$sel delete
mol delete top
exit



