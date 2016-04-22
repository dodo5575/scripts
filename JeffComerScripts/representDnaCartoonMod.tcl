# Author: Jeff Comer <jcomer2@illinois.edu>
proc representDna {segA segB colorList offset} {
    set backbone "C1' H1' C2' H2' H2'' C3' O3' H3' C4' O4' H4' C5' O5' H5' H5'' O1P O2P P"
    #set colorList {4 20 0 11 30}
  
    set n [molinfo top get numreps]
    for {set i 0} {$i < $n} {incr i} { mol delrep 0 top }

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

    set nc [llength $colorList]
    set offA [expr {($offset-$resA0)}]
    set offB [expr {(-$offset-$resB1)}]
    
    set count 0
    for {set i 0} {$i < $nc} {incr i} {
	set text "segname $segA and (resid+$offA+$i) % $nc == 0"
	
	mol representation NewCartoon 0.900000 20.000000 1.280000 0
	mol color ColorID [lindex $colorList $i]
	mol selection $text
	mol material Opaque
	mol addrep top
	incr count

	mol representation VDW 0.800000 15.000000
	mol color ColorID [lindex $colorList $i]
	mol selection "($text) and name C3'"
	mol material Opaque
	mol addrep top
	incr count
    }
 
    for {set i 0} {$i < $nc} {incr i} {
	set text "segname $segB and (resid+$offB-$i) % $nc == 0"
	
	mol representation NewCartoon 0.900000 20.000000 1.280000 0
	mol color ColorID [lindex $colorList $i]
	mol selection $text
	mol material Opaque
	mol addrep top
	incr count

	mol representation VDW 0.800000 15.000000
	mol color ColorID [lindex $colorList $i]
	mol selection "($text) and name C3'"
	mol material Opaque
	mol addrep top
	incr count
    }
}


color change rgb 0 0.2199999988079071 0.20999999344348907 0.8999999761581421
color change rgb 2 0.3499999940395355 0.3499999940395355 0.3499999940395355
color change rgb 3 0.36000001430511475 0.6299999952316284 0.25
color change rgb 4 0.8799999952316284 0.8799999952316284 0.30000001192092896
color change rgb 5 0.5 0.5 0.20000000298023224
color change rgb 6 0.7300000190734863 0.7300000190734863 0.7300000190734863
color change rgb 7 0.0 1.0 0.0
color change rgb 9 0.8999999761581421 0.5199999809265137 0.5199999809265137
color change rgb 11 0.6499999761581421 0.0 0.6499999761581421
color change rgb 12 0.5 0.8999999761581421 0.4000000059604645
color change rgb 13 0.8999999761581421 0.4000000059604645 0.699999988079071
color change rgb 14 0.5 0.30000001192092896 0.0
color change rgb 15 0.4699999988079071 0.4699999988079071 0.8999999761581421
color change rgb 17 0.8999999761581421 0.9399999976158142 0.800000011920929
color change rgb 18 0.550000011920929 0.8999999761581421 0.019999999552965164
color change rgb 19 0.0 0.8999999761581421 0.03999999910593033
color change rgb 20 0.0 0.8999999761581421 0.5
color change rgb 21 0.0 0.8799999952316284 1.0
color change rgb 22 0.0 0.6299999952316284 0.8899999856948853
color change rgb 23 0.019999999552965164 0.3799999952316284 0.6700000166893005
color change rgb 24 0.03999999910593033 0.10000000149011612 0.949999988079071
color change rgb 25 0.27000001072883606 0.0 0.9800000190734863
color change rgb 26 0.44999998807907104 0.0 0.8999999761581421
color change rgb 27 0.8999999761581421 0.6200000047683716 0.8999999761581421
color change rgb 28 0.6200000047683716 0.3100000023841858 1.0
color change rgb 29 0.9800000190734863 0.0 0.23000000417232513
color change rgb 30 0.8999999761581421 0.15000000596046448 0.15000000596046448
color change rgb 31 0.8899999856948853 0.3499999940395355 0.0
color change rgb 32 0.9599999785423279 0.7200000286102295 0.1599999964237213
