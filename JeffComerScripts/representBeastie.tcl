source createTubes.tcl

set n [molinfo top get numreps]
for {set i 0} {$i < $n} {incr i} { mol delrep 0 top }

mol representation Tube 1.3 40
mol color ColorID 0
mol material Opaque
mol selection "segname \"L.\""
mol addrep top

mol representation Tube 1.3 40
mol color ColorID 30
mol material Opaque
mol selection "segname \"R.\""
mol addrep top

#mol representation VDW 1.00000 21.000000
#mol color ColorID 0
#mol material Opaque
#mol selection "segname LA LB LC RA RB RC"
#mol addrep top

mol representation Tube 1.0 21 
mol color ColorID 2
mol material Opaque
mol selection "segname \"HR.*\""
mol addrep top

mol representation VDW 0.8 21
mol color ColorID 2
mol material Opaque
mol selection "segname \"HR.*\" and resid 0"
mol addrep top

color change rgb 0 0.2199999988079071 0.20999999344348907 0.8999999761581421
color change rgb 2 0.3499999940395355 0.3499999940395355 0.3499999940395355
color change rgb 3 0.36000001430511475 0.6299999952316284 0.25
color change rgb 4 0.8799999952316284 0.8799999952316284 0.30000001192092896
color change rgb 5 0.5 0.5 0.20000000298023224
color change rgb 6 0.7300000190734863 0.7300000190734863 0.7300000190734863
color change rgb 7 0.0 1.0 0.0
color change rgb 9 0.99 0.82 0.82
color change rgb 11 0.6499999761581421 0.0 0.6499999761581421
color change rgb 12 0.5 0.8999999761581421 0.4000000059604645
color change rgb 13 0.8999999761581421 0.4000000059604645 0.699999988079071
color change rgb 14 0.5 0.30000001192092896 0.0
color change rgb 15 0.82 0.82 0.99
color change rgb 17 0.8999999761581421 0.9399999976158142 0.800000011920929
color change rgb 18 0.550000011920929 0.8999999761581421 0.019999999552965164
color change rgb 19 0.1599999964237213 0.6600000262260437 0.27000001072883606
color change rgb 20 0.11999999731779099 0.6800000071525574 0.6700000166893005
color change rgb 21 0.0 0.8799999952316284 1.0
color change rgb 22 0.0 0.6299999952316284 0.8899999856948853
color change rgb 23 0.019999999552965164 0.3799999952316284 0.6700000166893005
color change rgb 24 0.03999999910593033 0.10000000149011612 0.949999988079071
color change rgb 25 0.27000001072883606 0.0 0.9800000190734863
color change rgb 26 0.44999998807907104 0.0 0.8999999761581421
color change rgb 27 0.8999999761581421 0.6200000047683716 0.8999999761581421
color change rgb 28 0.6200000047683716 0.3100000023841858 1.0
color change rgb 29 0.9800000190734863 0.0 0.23000000417232513
color change rgb 30 0.8 0.15 0.15
color change rgb 31 0.8899999856948853 0.3499999940395355 0.0
color change rgb 32 0.9599999785423279 0.7200000286102295 0.1599999964237213
