set all [atomselect top all]

$all moveby [vecinvert [measure center $all]]

set h0 "(segname T5 and resid  1 to 32) or (segname T6 and resid 33 to 64)"
set h1 "(segname T0 and resid 33 to 64) or (segname T1 and resid 33 to 64)"
set h2 "(segname T1 and resid  1 to 32) or (segname T3 and resid  1 to 32)"
set h3 "(segname T6 and resid  1 to 32) or (segname T4 and resid 33 to 64)"
set h4 "(segname T0 and resid  1 to 32) or (segname T2 and resid  1 to 32)"
set h5 "(segname T4 and resid  1 to 32) or (segname T5 and resid 33 to 64)"

set selText [list $h0 $h1 $h2 $h3 $h4 $h5]

foreach selT $selText {

    set sel [atomselect top $selT]
    $sel moveby [vecscale 0.1 [measure center $sel]]

    $sel delete
}


