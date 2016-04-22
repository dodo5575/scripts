# Fit the ends of molcules together.
# To use: vmd -dispdev text -e matchSegments.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set backboneAtoms {C1' H1' C2' H2' H2'' C3' H3' C4' H4' O4' C5' H5' H5''}
set selText "name $backboneAtoms"
set segListA {ADN0 ADN1 ADN2 ADNA ADN3 ADN4}
set segListB {BDN0 BDN1 BDN2 BDNA BDN3 BDN4}
set resA0 1
set resA1 58
set resB0 111
set resB1 54
# Input:
set psf double_much.psf
set pdb double_much.pdb
# Output:
set outPrefix double_long

mol load psf $psf pdb $pdb

set n [llength $segListA]

for {set i 1} {$i < $n} {incr i} {
    set segA0 [lindex $segListA [expr $i-1]]
    set segB0 [lindex $segListB [expr $i-1]]
    set sel0 [atomselect top "($selText) and ((segname $segA0 and resid $resA1) or (segname $segB0 and resid $resB1))"]
    
    set segA1 [lindex $segListA $i]
    set segB1 [lindex $segListB $i]
    set sel1 [atomselect top "($selText) and ((segname $segA1 and resid $resA0) or (segname $segB1 and resid $resB0))"]
    set all1 [atomselect top "segname $segA1 $segB1"]

    
    set m [measure fit $sel1 $sel0]
    $all1 move $m

    $sel0 delete
    $sel1 delete
    $all1 delete
}

set all [atomselect top all]
$all writepdb $outPrefix.pdb
$all writepsf $outPrefix.psf
$all delete

mol delete top
exit



