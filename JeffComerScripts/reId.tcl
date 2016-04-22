set psf fit1.psf
set pdb fit1.pdb
set outPdb fit2.pdb

proc reId {seg res0 res1 newRes0} {
    set selList {}
    for {set i $res0} {$i <= $res1} {incr i} {
	lappend selList [atomselect top "segname $seg and resid $i"]
    }

    foreach s $selList {
	set res [lindex [$s get resid] 0]
	$s set resid [expr $res-$res0+$newRes0]
	$s delete
    }
    return
}

mol load psf $psf pdb $pdb
reId ADNA 2 18 102
reId DNA2 4 24 119
reId DNA1 1 23 [expr 101-22]

set sel [atomselect top "(segname ADNA and resid 102 to 118) or (segname BDNA and resid 99 to 110) or (segname DNA1 and resid [expr 101-22] to 101) or (segname DNA2 and resid 119 to 139)"]
$sel writepdb $outPdb



