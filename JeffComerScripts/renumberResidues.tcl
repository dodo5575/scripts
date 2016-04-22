# Author: Jeff Comer <jcomer2@illinois.edu>
set startResidue 1
# Input:
set inName ss_triplet_conform3a
# Output:
set outName ss_triplet_residue

mol load psf $inName.psf pdb $inName.pdb
set all [atomselect top all]
set segList [lsort -unique [$all get segname]]

set count 0
foreach seg $segList {
    set sel [atomselect top "segname $seg"]
    set resList [lsort -unique -integer [$sel get resid]]
    $sel delete

    # Select the residues.
    set selList {}
    foreach res $resList {
	lappend selList [atomselect top "segname $seg and resid $res"]
    } 

    set r $startResidue
    foreach sel $selList res $resList {
	$sel set resid $r
	puts "$seg:$res -> $seg:$r"
	incr r
	incr count
    }
}

puts "Renumbered $count residues."

$all writepsf $outName.psf
$all writepdb $outName.pdb
exit
