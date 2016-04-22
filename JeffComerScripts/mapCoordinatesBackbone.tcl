# This script combines two pdb/psf pairs.
# Use with: vmd -dispdev text -e combine.tcl
# Author: Jeff Comer <jcomer2@uiuc.edu>

set backboneAtoms "C1' H1' C2' H2' H2'' C3' O3' H3' C4' O4' H4' C5' O5' H5' H5'' O1P O2P P"
set selText "segname ADNA BDNA and name $backboneAtoms"
# Input:
set coordPdb cut_dna.pdb
set mapPsf dna_A-T_fix.psf
set mapPdb dna_A-T_fix.pdb
# Output:
set finalPdb dna_A-T_map.pdb

#set coordMol [mol load psf $coordPsf pdb $coordPdb]
set coordMol [mol load pdb $coordPdb]
set mapMol [mol load psf $mapPsf pdb $mapPdb]

set mapSel [atomselect $mapMol $selText]

foreach zero {0} {
    set nameList [$mapSel get name]
    set resList [$mapSel get resid]
    set segList [$mapSel get segname]
}
$mapSel delete

set res0 0
foreach nam $nameList res $resList seg $segList {
    if {$res0 != $res} {
	puts "$seg:$res"
	set res0 $res
    }

    set mSel [atomselect $mapMol "segname $seg and resid $res and name $nam"]
    set cSel [atomselect $coordMol "segname $seg and resid $res and name $nam"]

    if {[$cSel num] == 0} {
	puts "WARNING! Cannot find $seg:$res:$nam in $coordPdb."
    } else {
        $mSel set {x y z} [$cSel get {x y z}]
    }
    

    $mSel delete
    $cSel delete
}

# Save the mapped molecule.
set all [atomselect $mapMol all]
$all writepdb $finalPdb


mol delete top
mol delete top
exit


