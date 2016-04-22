# Use with: vmd -dispdev text -e integrateSegments.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set segList {ADNA ADNB ADNC ADND}
#Input:
set psf cont_dna_KCl_types.psf
set pdb cont_dna_KCl.pdb
#Output:
set outName cont_dna_int
set finalPsf ${outName}.psf
set finalPdb ${outName}.pdb

set resId 1
set segName [lindex $segList 0]
mol load psf $psf pdb $pdb

foreach seg $segList {
    puts "Integrating segment $seg into segment $segName..."
    set sel [atomselect top "segname $seg"]
    set resList [lsort -unique -integer [$sel get resid]]
    $sel delete

    # Renumber the residues for each segment.
    foreach res $resList {
	set s [atomselect top "segname $seg and resid $res"]
	$s set resid $resId
	$s set segname $segName
	$s delete

	incr resId
    }
}

set all [atomselect top all]
$all writepsf $finalPsf
$all writepdb $finalPdb

mol delete top
exit


