# This script combines two pdb/psf pairs.
# Use with: vmd -dispdev text -e combine.tcl
# Author: Jeff Comer <jcomer2@uiuc.edu>

set backboneAtoms "C1' H1' C2' H2' H2'' C3' O3' H3' C4' O4' H4' C5' O5' H5' H5'' O1P O2P P"
set selText "nucleic"
set resShift 7
# Input:
set coordPdb total_frame64.pdb
set mapPsf ln_polyA_frame5.psf
set mapPdb total_frame64.pdb
# Output:
set finalPdb shift_total_frame64.pdb

set coordMol [mol load pdb $coordPdb]
set mapMol [mol load psf $mapPsf pdb $mapPdb]

set mapSel [atomselect $mapMol $selText]

foreach zero {0} {
    set resList [lsort -unique -integer [$mapSel get residue]]
    set atomList [lsort -integer -index 1 [$mapSel get {segname resid name residue}]]
}
$mapSel delete

set res0 0
foreach atom $atomList {
    foreach {seg res nam residue} $atom { break }

    # Find the shifted the residue.
    set mapResidueInd [lsearch $resList $residue]
    set coordResidueInd [expr {$mapResidueInd + $resShift}]
    if {$coordResidueInd >= [llength $resList]} {
	set coordResidueInd [expr {$coordResidueInd - [llength $resList]}]
    }
    set coordResidue [lindex $resList $coordResidueInd]

    if {$res0 != $res} {
	puts "$mapResidueInd"
	puts "$seg:$res residue $residue -> $coordResidue"
	set res0 $res
    }

    set mSel [atomselect $mapMol "residue $residue and name $nam"]
    set cSel [atomselect $coordMol "residue $coordResidue and name $nam"]

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
