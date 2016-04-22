# Use with: vmd -dispdev text -e changeIds.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set psf clamp_dna.psf
set pdb clamp_dna.pdb
set outName clamp_dna_id

set seg BDNA

mol load psf $psf pdb $pdb
set sel [atomselect top "segname $seg"]

set resList [lsort -unique -integer [$sel get resid]]

set selList {}
foreach res $resList {
    lappend selList [atomselect top "segname $seg and resid $res"]
}

foreach res $resList s $selList {
    set res1 [expr $res + 20]
    if {$res1 > 30} {
	set res1 [expr $res1 - 30]
    }

    $s set resid $res1
    puts "Changed residue $res to $res1."
}

set all [atomselect top all]
$all writepsf $outName.psf
$all writepdb $outName.pdb
$all delete

foreach s $selList {
    $s delete
}
mol delete top
exit
