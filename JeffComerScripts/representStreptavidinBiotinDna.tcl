set backbone "C1' H1' C2' H2' H2'' C3' O3' H3' C4' O4' H4' C5' O5' H5' H5'' O1P O2P P"
#mol delrep all top
set n [molinfo top get numreps]
for {set i 0} {$i < $n} {incr i} { mol delrep 0 top }

set segList1 {STPA STPB STPC STPD}
set segList {ADNA BDNA}
set darkList {30 0}
set lightList {9 15}

foreach seg $segList dark $darkList light $lightList {
    mol representation NewCartoon 0.900000 20.000000 1.280000 0
    mol color ColorID $dark
    mol selection "segname $seg"
    mol material Opaque
    mol addrep top

    mol representation VDW 0.900000 15.000000
    mol color ColorID $dark
    mol selection "segname $seg and name C3'"
    mol material Opaque
    mol addrep top
}

mol representation VDW 0.900000 15.000000
mol color ColorID 4
mol selection "segname [lrange $segList1 0 1]"
mol material Opaque
mol addrep top

mol representation VDW 0.900000 15.000000
mol color ColorID 19
mol selection "segname [lrange $segList1 2 3]"
mol material Opaque
mol addrep top

mol representation VDW 0.900000 15.000000
mol color ColorID 11
mol selection "segname BTN"
mol material Opaque
mol addrep top

mol representation VDW 0.900000 15.000000
mol color ColorID 30
mol selection "segname ADNA and resid 49 and name P3 O1P3 O2P3 O3T O3' C3'"
mol material Opaque
mol addrep top

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
color change rgb 19 0.16 0.66 0.27
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
