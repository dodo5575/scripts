set backbone "C1' H1' C2' H2' H2'' C3' O3' H3' C4' O4' H4' C5' O5' H5' H5'' O1P O2P P"

set selA [atomselect top "segname $segA"]
set resListA [lsort -unique -integer [$selA get resid]]
$selA delete
set resA0 [lindex $resListA 0]
set resA1 [lindex $resListA end]
set numA [expr {$resA1 - $resA0}]

set selB [atomselect top "segname $segB"]
set resListB [lsort -unique -integer [$selB get resid]]
$selB delete
set resB0 [lindex $resListB 0]
set resB1 [lindex $resListB end]
set numB [expr {$resB1 - $resB0}]

set nBasepairs $numA
if {$numB < $numA} { set nBasepairs $numB }
