# This script combines two pdb/psf pairs.
# Use with: vmd -dispdev text -e combine.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set selText "(segname ADNA and not resid 28) or (segname BDNA and not resid 35)"
# Input:
set coordPsf bam_all.psf
set coordPdb bam_all.pdb
set mapPsf dna_bam_spec.psf
set mapPdb dna_bam_spec.pdb
# Output:
set finalPdb dna_bam_spec_mapped.pdb

set coordStart 1
set mapStart 1
set nResidues 62

set coordMol [mol load psf $coordPsf pdb $coordPdb]
set mapMol [mol load psf $mapPsf pdb $mapPdb]

for {set i 0} {$i < $nResidues} {incr i} {
    set cx [expr $coordStart + $i]
    set coordSel [atomselect $coordMol "($selText) and resid $cx"]
    set coordResName [lindex [$coordSel get resname] 0]
 
    set mx [expr $mapStart + $i]
    set mapSel [atomselect $mapMol "($selText) and resid $mx"]
    set mapResName [lindex [$mapSel get resname] 0]

    if {![string equal $coordResName $mapResName]} {
	puts "ERROR: Residues ${coordResName}${cx} and ${mapResName}${mx} do not agree!"
	exit
    }

    set nameList [$coordSel get name]
    set nameList [$mapSel get name]
    $coordSel delete
    $mapSel delete
    
    foreach name $nameList {
	set cs [atomselect $coordMol "($selText) and resid $cx and name $name"]
	set cr [$cs get {x y z}]
	$cs delete
	
	set ms [atomselect $mapMol "($selText) and resid $mx and name $name"]
	if {[$ms num] == 0} {
	    puts "ERROR: Atom name $name is not found in map molecule."
	    exit
	}
	$ms set {x y z} $cr
	$ms delete
    }
}

# Save the mapped molecule.
set all [atomselect $mapMol all]
$all writepdb $finalPdb


mol delete top
mol delete top
exit



