set psf thin_dna_KCl.psf
set coor thin_dna_KCl.pdb
set outPrefix thin_dna_resid

proc reId1 {seg uniName} {
    set sel [atomselect top "segname $seg and name $uniName"]
    set resList [lsort -integer -index 1 [$sel get {resid index}]]
    $sel delete

    set selList {}
    foreach item $resList {
	lappend selList [atomselect top "segname $seg and resid [lindex $item 0]"]
    }

    set res 1
    foreach s $selList item $resList {
 	$s set resid $res
	puts "Changed residue [lindex $item 0] to $res ([$s num] atoms)."
	incr res
	$s delete
    }

    return
}

mol load psf $psf
mol addfile $coor
set all [atomselect top all]
reId1 BDNA P
$all writepsf $outPrefix.psf
$all writepdb $outPrefix.pdb
exit
