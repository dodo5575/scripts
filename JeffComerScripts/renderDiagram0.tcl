# Author: Jeff Comer <jcomer2@illinois.edu>

source $env(HOME)/scripts/vector.tcl
set pi [expr {4.0*atan(1.0)}]

# Construct a pdb line.
proc makePdbLineFull {index segName resId name resName r} {
    set template "ATOM     42  HN  ASN X   4     -41.083  17.391  50.684  0.00  0.00      P1    "

    foreach {x y z} $r {break}
    set record "ATOM  "
    set si [string range [format "     %5i " $index] end-5 end]
    if {[string length $name] < 4} {
	set name [string range " $name    " 0 3]
    } else {
	set name [string range $name 0 3]
    }
    set resName [string range " $resName    " 0 3]
    set temp0 " [string index $segName 0]"
    set resId [string range "    $resId"  end-3 end]
    set temp1 [string range $template  26 29]
    set sx [string range [format "       %8.3f" $x] end-7 end]
    set sy [string range [format "       %8.3f" $y] end-7 end]
    set sz [string range [format "       %8.3f" $z] end-7 end]
    set temp2 [string range $template 54 71]
    set segName [string range "$segName    "  0 3]
    set tempEnd [string range $template 76 end]

    # Construct the pdb line.
    return "${record}${si}${name}${resName}${temp0}${resId}${temp1}${sx}${sy}${sz}${temp2}${segName}${tempEnd}"
}


proc makePdbHeader {{comment "generated pdb"}} {
    # Write the pdb header.
    return "REMARK $comment"
}

proc makeTubePdb {outPdb posListList} {
    set out [open $outPdb w]
    puts $out [makePdbHeader "tubes"]

    set res 0
    foreach posList $posListList {
	set i 0
	foreach r $posList {
	    puts $out [makePdbLineFull $i TU${res} $i CA ALA $r]
	    incr i
	}
	incr res
    }
    close $out
}

proc makeComboTubePdb {outPdb posListList segList} {
    set out [open $outPdb w]
    puts $out [makePdbHeader "tubes"]

    foreach posList $posListList seg $segList {
	set i 0
	foreach r $posList {
	    puts $out [makePdbLineFull $i $seg $i CA ALA $r]
	    incr i
	}
    }
    close $out
}

proc interpDirection {r0 r1 factor} {
    set u0 [vecUnit $r0]
    set u1 [vecUnit $r1]
    set n [vecCross $u0 $u1]
    set phi [expr {acos([vecDot $u0 $u1])}]
    set t [expr {$phi*$factor}]
    
    set ct [expr {cos($t)}]
    set st [expr {sin($t)}]
    set u [vecAdd [vecScale $ct $u0] [vecScale $st $u1]]
    return $u
}

proc rotVectorPerp {r axis theta} {
    set ct [expr {cos($theta)}]
    set st [expr {sin($theta)}]
    set u [vecUnit $r]
    set n [vecCross $a $u]
    return [vecAdd [vecScale $ct $r] [vecScale $st $n]]
}

proc randomPerp {unit} {
    global pi
    set n [vecCross $unit [list 0 0 1]]
    if {[vecLength $n] == 0.0} { set n [vecCross $unit [list 1 0 0]] }
    set n [vecUnit $n]
    set b [vecCross $unit $n]

    set theta [expr {2.0*$pi*rand()}]
    set ct [expr {cos($theta)}]
    set st [expr {sin($theta)}]
    return [vecAdd [vecScale $ct $n] [vecScale $st $b]]
}

proc upPerp {unit up} {
    global pi
    set proj [vecScale [vecDot $up $unit] $unit]
    set n [vecSub $up $proj]
    if {[vecLength $n] == 0.0} { set n [vecCross $unit [list 0.1 0.2 0.3]] }
    set n [vecUnit $n]
    return $n
}

proc makeStrand {num initPos initDir len angleMax biasVector biasFactor zMin} {
    set r0 $initPos
    set d0 $initDir
    set u0 [vecUnit $d0]
    #set n0 [vecCross $d0 [list 0 0 1]]
    #if {[vecLength $n0] == 0.0} { set n0 [vecCross $d0 [list 1 0 0]] }
    set n0 [randomPerp $u0]
    
    set pi [expr {4.0*atan(1.0)}]
    set phiMax [expr {$angleMax/180.0*$pi}]

    # Check for infinite loop.
    if {[lindex $initPos 2] < $zMin} { return {} }

    set posList [list $r0]
    set j 0
    while {$j < $num} {
	# Pick a random angle of the direction as large as phiMax.
	set phi [expr {$phiMax*2.0*(rand()-0.5)}]
	# Pick the azimuthal direction at random.
	set theta [expr {2.0*$pi*rand()}]
	
	# Compute the new chain direction.
	set ct [expr {cos($phi)}]
	set st [expr {sin($phi)}]
	set u [vecAdd [vecScale $ct $u0] [vecScale $st $n0]]

	# Add the bias.
	set u [interpDirection $u $biasVector $biasFactor]
	set r [vecAdd $r0 [vecScale $len $u]]
	
	if {[lindex $r 2] > $zMin} {
	    lappend posList $r
	    set r0 $r
	    set u0 $u
	    set n0 [randomPerp $u0]
	    incr j
	}
    }
    return $posList
}

proc makeChain {pos0 pos1 lenFactor num up} {
    # We need an odd number.
    if {$num % 2 == 0} {set num [expr {$num + 1}]}
    if {$lenFactor < 1.01} { set lenFactor 1.01 }

    # We need to make a curve of length $len that goes through $pos0 and $pos1.
    set d [vecSub $pos1 $pos0]
    set dMag [vecLength $d]
    set len [expr {$lenFactor*$dMag}]
    set tang [vecUnit $d]
    set norm [upPerp $tang $up]

    set s [expr {0.5*$len}]
    set l [expr {0.5*[vecLength $d]}]
    set height [expr {sqrt($s*$s - $l*$l)}]
    set mid [expr {($num-1)/2}]
    set dy [expr {$height/$mid}]
    set dx [expr {$l/$mid}]

    # Start the position list.
    set posList [list $pos0]

    # Make the first leg of the triangle.
    for {set i 1} {$i <= $mid} {incr i} {
	set vx [vecScale [expr {$dx*$i}] $tang]
	set vy [vecScale [expr {$dy*$i}] $norm]
	set r [vecAdd $pos0 [vecAdd $vx $vy]]
	lappend posList $r
    }

    # Make the second leg of the triangle.
    for {set i [expr {$mid+1}]} {$i < [expr {$num-1}]} {incr i} {
	set j [expr {2*$mid-$i}]
	set vx [vecScale [expr {$dx*$i}] $tang]
	set vy [vecScale [expr {$dy*$j}] $norm]
	set r [vecAdd $pos0 [vecAdd $vx $vy]]
	lappend posList $r
    }

    # Add the last point.
    lappend posList $pos1
    return $posList
}

proc monteCarloSteps {point n } {
    for {set i 0} {$i < $n} {incr i} {
	set point [monteCarloRotation $point]
    }
    return $point
}

proc monteCarloGood {point n zMin} {
    for {set i 0} {$i < $n} {incr i} {
	set p [monteCarloRotation $point]
	
	# Check for collisions.
	set good 1
	foreach r $p {
	    if {[lindex $r 2] < $zMin} {
		set good 0
		break
	    }
	}
	if {$good} { set point $p }
    }
    return $point
}

proc monteCarloRotation {posList} {
    global pi
    set n [llength $posList]

    # Pick two points about which to perform the rotation.
    set pick0 [expr {int(floor(rand()*$n))}]
    set pick1 [expr {int(floor(rand()*$n))}]
    if {$pick0 == $pick1} { return $posList }
    if {$pick1 < $pick0} { 
	set tmp $pick0
	set pick0 $pick1
	set pick1 $tmp
    }

    # Get the angle and axis.
    set r0 [lindex $posList $pick0]
    set r1 [lindex $posList $pick1]
    set d [vecSub $r1 $r0]
    set angle [expr {2.0*$pi*rand()}]
    set ct [expr {cos($angle)}]
    set st [expr {sin($angle)}]
    set axis [vecUnit [vecSub $r1 $r0]]

    # Rotate each point.
    for {set i [expr {$pick0+1}]} {$i < $pick1} {incr i} {
	set r [vecSub [lindex $posList $i] $r0]
	set z [vecScale [vecDot $r $axis] $axis]
	set y [vecCross $axis $r]

	set dx [vecScale $ct [vecSub $r $z]]
	set dy [vecScale $st $y]
	set dz $z
	lset posList $i [vecAdd [vecAdd [vecAdd $dx $dy] $dz] $r0]
    }
    return $posList
}

proc drawChargeBall {center rotation radius sign baseColor symbolColor} {
    global pi
    set meridians 19
    set parallels 15

    set greenwich [expr {$meridians/2}]
    set equator [expr {$parallels/2+1}]

    set dTheta [expr {2.0*$pi/$meridians}]
    set dPhi [expr {$pi/$parallels}]

    # Fill in the patches on the sphere.
    for {set p 2} {$p < $parallels} {incr p} {
	set phi0 [expr {($p-1)*$dPhi}]
	set phi1 [expr {$p*$dPhi}]
	for {set m 0} {$m < $meridians} {incr m} {
	    set theta0 [expr {$m*$dTheta}]
	    set theta1 [expr {($m+1)*$dTheta}]
	    
	    set n(nw) [spherePos $phi0 $theta0]
	    set n(ne) [spherePos $phi0 $theta1]
	    set n(sw) [spherePos $phi1 $theta0]
	    set n(se) [spherePos $phi1 $theta1]
	    
	    foreach d {nw ne sw se} {
		set r($d) [vecScale $radius $n($d)]
	    }

	    # Do the transformation.
	    foreach d {nw ne sw se} {
		set n($d) [vecTransform $rotation $n($d)]
		set r($d) [vecAdd [vecTransform $rotation $r($d)] $center]
	    }

	    # Add the horizontal bar.
	    if {abs($greenwich-$m) < 3 && $equator == $p} {
		graphics top color $symbolColor
	    } else {
		graphics top color $baseColor
	    }
	    
	    # For a positive charge, add the vertical bar.
	    if {$sign > 0.0} {
		if {abs($equator-$p) < 4 && $greenwich == $m} {
		    graphics top color $symbolColor
		}
	    }
	    
	    graphics top trinorm $r(nw) $r(se) $r(ne) $n(nw) $n(se) $n(ne)
	    graphics top trinorm $r(nw) $r(sw) $r(se) $n(nw) $n(sw) $n(se)
	}
    }

    # Do the poles.
    graphics top color $baseColor
    set phiA $dPhi
    set phiB [expr {($parallels-1)*$dPhi}]
    
    set n(a) {0.0 0.0 1.0}
    set n(b) {0.0 0.0 -1.0}
    set r(a) [vecScale $radius $n(a)]
    set r(b) [vecScale $radius $n(b)]	
    # Do the transformation.
    foreach d {a b} {
	set n($d) [vecTransform $rotation $n($d)]
	set r($d) [vecAdd [vecTransform $rotation $r($d)] $center]
    }

    for {set m 0} {$m < $meridians} {incr m} {
	set theta0 [expr {$m*$dTheta}]
	set theta1 [expr {($m+1)*$dTheta}]
	
	set n(nw) [spherePos $phiA $theta0]
	set n(ne) [spherePos $phiA $theta1]
	
	set n(sw) [spherePos $phiB $theta0]
	set n(se) [spherePos $phiB $theta1]
	
	foreach d {nw ne sw se} {
	    set r($d) [vecScale $radius $n($d)]
	}

	# Do the transformation.
	foreach d {nw ne sw se} {
	    set n($d) [vecTransform $rotation $n($d)]
	    set r($d) [vecAdd [vecTransform $rotation $r($d)] $center]
	}

	graphics top trinorm $r(a) $r(nw) $r(ne) $n(a) $n(nw) $n(ne)
	graphics top trinorm $r(b) $r(se) $r(sw) $n(b) $n(se) $n(sw)
    }
    return
}

proc drawChargeBall {center rotation radius sign baseColor symbolColor} {
    global pi
    set meridians 19
    set parallels 15

    set greenwich [expr {$meridians/2}]
    set equator [expr {$parallels/2+1}]

    set dTheta [expr {2.0*$pi/$meridians}]
    set dPhi [expr {$pi/$parallels}]

    # Fill in the patches on the sphere.
    for {set p 2} {$p < $parallels} {incr p} {
	set phi0 [expr {($p-1)*$dPhi}]
	set phi1 [expr {$p*$dPhi}]
	for {set m 0} {$m < $meridians} {incr m} {
	    set theta0 [expr {$m*$dTheta}]
	    set theta1 [expr {($m+1)*$dTheta}]
	    
	    set n(nw) [spherePos $phi0 $theta0]
	    set n(ne) [spherePos $phi0 $theta1]
	    set n(sw) [spherePos $phi1 $theta0]
	    set n(se) [spherePos $phi1 $theta1]
	    
	    foreach d {nw ne sw se} {
		set r($d) [vecScale $radius $n($d)]
	    }

	    # Do the transformation.
	    foreach d {nw ne sw se} {
		set n($d) [vecTransform $rotation $n($d)]
		set r($d) [vecAdd [vecTransform $rotation $r($d)] $center]
	    }

	    # Add the horizontal bar.
	    if {abs($greenwich-$m) < 3 && $equator == $p} {
		graphics top color $symbolColor
	    } else {
		graphics top color $baseColor
	    }
	    
	    # For a positive charge, add the vertical bar.
	    if {$sign > 0.0} {
		if {abs($equator-$p) < 4 && $greenwich == $m} {
		    graphics top color $symbolColor
		}
	    }
	    
	    graphics top trinorm $r(nw) $r(se) $r(ne) $n(nw) $n(se) $n(ne)
	    graphics top trinorm $r(nw) $r(sw) $r(se) $n(nw) $n(sw) $n(se)
	}
    }

    # Do the poles.
    graphics top color $baseColor
    set phiA $dPhi
    set phiB [expr {($parallels-1)*$dPhi}]
    
    set n(a) {0.0 0.0 1.0}
    set n(b) {0.0 0.0 -1.0}
    set r(a) [vecScale $radius $n(a)]
    set r(b) [vecScale $radius $n(b)]	
    # Do the transformation.
    foreach d {a b} {
	set n($d) [vecTransform $rotation $n($d)]
	set r($d) [vecAdd [vecTransform $rotation $r($d)] $center]
    }

    for {set m 0} {$m < $meridians} {incr m} {
	set theta0 [expr {$m*$dTheta}]
	set theta1 [expr {($m+1)*$dTheta}]
	
	set n(nw) [spherePos $phiA $theta0]
	set n(ne) [spherePos $phiA $theta1]
	
	set n(sw) [spherePos $phiB $theta0]
	set n(se) [spherePos $phiB $theta1]
	
	foreach d {nw ne sw se} {
	    set r($d) [vecScale $radius $n($d)]
	}

	# Do the transformation.
	foreach d {nw ne sw se} {
	    set n($d) [vecTransform $rotation $n($d)]
	    set r($d) [vecAdd [vecTransform $rotation $r($d)] $center]
	}

	graphics top trinorm $r(a) $r(nw) $r(ne) $n(a) $n(nw) $n(ne)
	graphics top trinorm $r(b) $r(se) $r(sw) $n(b) $n(se) $n(sw)
    }
    return
}

proc drawChargeBall {center rotation radius sign baseColor symbolColor} {
    global pi
    set meridians 19
    set parallels 15

    set greenwich [expr {$meridians/2}]
    set equator [expr {$parallels/2+1}]

    set dTheta [expr {2.0*$pi/$meridians}]
    set dPhi [expr {$pi/$parallels}]

    # Fill in the patches on the sphere.
    for {set p 2} {$p < $parallels} {incr p} {
	set phi0 [expr {($p-1)*$dPhi}]
	set phi1 [expr {$p*$dPhi}]
	for {set m 0} {$m < $meridians} {incr m} {
	    set theta0 [expr {$m*$dTheta}]
	    set theta1 [expr {($m+1)*$dTheta}]
	    
	    set n(nw) [spherePos $phi0 $theta0]
	    set n(ne) [spherePos $phi0 $theta1]
	    set n(sw) [spherePos $phi1 $theta0]
	    set n(se) [spherePos $phi1 $theta1]
	    
	    foreach d {nw ne sw se} {
		set r($d) [vecScale $radius $n($d)]
	    }

	    # Do the transformation.
	    foreach d {nw ne sw se} {
		set n($d) [vecTransform $rotation $n($d)]
		set r($d) [vecAdd [vecTransform $rotation $r($d)] $center]
	    }

	    # Add the horizontal bar.
	    if {abs($greenwich-$m) < 3 && $equator == $p} {
		graphics top color $symbolColor
	    } else {
		graphics top color $baseColor
	    }
	    
	    # For a positive charge, add the vertical bar.
	    if {$sign > 0.0} {
		if {abs($equator-$p) < 4 && $greenwich == $m} {
		    graphics top color $symbolColor
		}
	    }
	    
	    graphics top trinorm $r(nw) $r(se) $r(ne) $n(nw) $n(se) $n(ne)
	    graphics top trinorm $r(nw) $r(sw) $r(se) $n(nw) $n(sw) $n(se)
	}
    }

    # Do the poles.
    graphics top color $baseColor
    set phiA $dPhi
    set phiB [expr {($parallels-1)*$dPhi}]
    
    set n(a) {0.0 0.0 1.0}
    set n(b) {0.0 0.0 -1.0}
    set r(a) [vecScale $radius $n(a)]
    set r(b) [vecScale $radius $n(b)]	
    # Do the transformation.
    foreach d {a b} {
	set n($d) [vecTransform $rotation $n($d)]
	set r($d) [vecAdd [vecTransform $rotation $r($d)] $center]
    }

    for {set m 0} {$m < $meridians} {incr m} {
	set theta0 [expr {$m*$dTheta}]
	set theta1 [expr {($m+1)*$dTheta}]
	
	set n(nw) [spherePos $phiA $theta0]
	set n(ne) [spherePos $phiA $theta1]
	
	set n(sw) [spherePos $phiB $theta0]
	set n(se) [spherePos $phiB $theta1]
	
	foreach d {nw ne sw se} {
	    set r($d) [vecScale $radius $n($d)]
	}

	# Do the transformation.
	foreach d {nw ne sw se} {
	    set n($d) [vecTransform $rotation $n($d)]
	    set r($d) [vecAdd [vecTransform $rotation $r($d)] $center]
	}

	graphics top trinorm $r(a) $r(nw) $r(ne) $n(a) $n(nw) $n(ne)
	graphics top trinorm $r(b) $r(se) $r(sw) $n(b) $n(se) $n(sw)
    }
    return
}

proc drawChargeBall {center rotation radius sign baseColor symbolColor} {
    global pi
    set meridians 19
    set parallels 15

    set greenwich [expr {$meridians/2}]
    set equator [expr {$parallels/2+1}]

    set dTheta [expr {2.0*$pi/$meridians}]
    set dPhi [expr {$pi/$parallels}]

    # Fill in the patches on the sphere.
    for {set p 2} {$p < $parallels} {incr p} {
	set phi0 [expr {($p-1)*$dPhi}]
	set phi1 [expr {$p*$dPhi}]
	for {set m 0} {$m < $meridians} {incr m} {
	    set theta0 [expr {$m*$dTheta}]
	    set theta1 [expr {($m+1)*$dTheta}]
	    
	    set n(nw) [spherePos $phi0 $theta0]
	    set n(ne) [spherePos $phi0 $theta1]
	    set n(sw) [spherePos $phi1 $theta0]
	    set n(se) [spherePos $phi1 $theta1]
	    
	    foreach d {nw ne sw se} {
		set r($d) [vecScale $radius $n($d)]
	    }

	    # Do the transformation.
	    foreach d {nw ne sw se} {
		set n($d) [vecTransform $rotation $n($d)]
		set r($d) [vecAdd [vecTransform $rotation $r($d)] $center]
	    }

	    # Add the horizontal bar.
	    if {abs($greenwich-$m) < 3 && $equator == $p} {
		graphics top color $symbolColor
	    } else {
		graphics top color $baseColor
	    }
	    
	    # For a positive charge, add the vertical bar.
	    if {$sign > 0.0} {
		if {abs($equator-$p) < 4 && $greenwich == $m} {
		    graphics top color $symbolColor
		}
	    }
	    
	    graphics top trinorm $r(nw) $r(se) $r(ne) $n(nw) $n(se) $n(ne)
	    graphics top trinorm $r(nw) $r(sw) $r(se) $n(nw) $n(sw) $n(se)
	}
    }

    # Do the poles.
    graphics top color $baseColor
    set phiA $dPhi
    set phiB [expr {($parallels-1)*$dPhi}]
    
    set n(a) {0.0 0.0 1.0}
    set n(b) {0.0 0.0 -1.0}
    set r(a) [vecScale $radius $n(a)]
    set r(b) [vecScale $radius $n(b)]	
    # Do the transformation.
    foreach d {a b} {
	set n($d) [vecTransform $rotation $n($d)]
	set r($d) [vecAdd [vecTransform $rotation $r($d)] $center]
    }

    for {set m 0} {$m < $meridians} {incr m} {
	set theta0 [expr {$m*$dTheta}]
	set theta1 [expr {($m+1)*$dTheta}]
	
	set n(nw) [spherePos $phiA $theta0]
	set n(ne) [spherePos $phiA $theta1]
	
	set n(sw) [spherePos $phiB $theta0]
	set n(se) [spherePos $phiB $theta1]
	
	foreach d {nw ne sw se} {
	    set r($d) [vecScale $radius $n($d)]
	}

	# Do the transformation.
	foreach d {nw ne sw se} {
	    set n($d) [vecTransform $rotation $n($d)]
	    set r($d) [vecAdd [vecTransform $rotation $r($d)] $center]
	}

	graphics top trinorm $r(a) $r(nw) $r(ne) $n(a) $n(nw) $n(ne)
	graphics top trinorm $r(b) $r(se) $r(sw) $n(b) $n(se) $n(sw)
    }
    return
}

proc drawChargeBallClear {center rotation radius sign symbolColor} {
    global pi
    set meridians 19
    set parallels 15

    set greenwich [expr {$meridians/2}]
    set equator [expr {$parallels/2+1}]

    set dTheta [expr {2.0*$pi/$meridians}]
    set dPhi [expr {$pi/$parallels}]

    # Fill in the patches on the sphere.
    for {set p 2} {$p < $parallels} {incr p} {
	set phi0 [expr {($p-1)*$dPhi}]
	set phi1 [expr {$p*$dPhi}]
	for {set m 0} {$m < $meridians} {incr m} {
	    set theta0 [expr {$m*$dTheta}]
	    set theta1 [expr {($m+1)*$dTheta}]
	    
	    set n(nw) [spherePos $phi0 $theta0]
	    set n(ne) [spherePos $phi0 $theta1]
	    set n(sw) [spherePos $phi1 $theta0]
	    set n(se) [spherePos $phi1 $theta1]
	    
	    foreach d {nw ne sw se} {
		set r($d) [vecScale $radius $n($d)]
	    }

	    # Do the transformation.
	    foreach d {nw ne sw se} {
		set n($d) [vecTransform $rotation $n($d)]
		set r($d) [vecAdd [vecTransform $rotation $r($d)] $center]
	    }

	    # Draw the symbol 
	    graphics top color $symbolColor
	    if {$sign > 0.0} {
		if {(abs($greenwich-$m) < 3 && $equator == $p) 
		    || (abs($equator-$p) < 4 && $greenwich == $m)} {
		    graphics top trinorm $r(nw) $r(se) $r(ne) $n(nw) $n(se) $n(ne)
		    graphics top trinorm $r(nw) $r(sw) $r(se) $n(nw) $n(sw) $n(se)
		} 
	    } else {
		if {(abs($greenwich-$m) < 3 && $equator == $p)} {
		    graphics top trinorm $r(nw) $r(se) $r(ne) $n(nw) $n(se) $n(ne)
		    graphics top trinorm $r(nw) $r(sw) $r(se) $n(nw) $n(sw) $n(se)
		} 
	    }
	}
    }
    return
}

proc drawPill {center len diam {mole top}} {
    set l [expr {0.5*$len}]
    set s [expr {0.5*$diam}]

    set r0 [vecSub $center [list $l 0.0 0.0]]
    set r1 [vecAdd $center [list $l 0.0 0.0]]

    graphics $mole cylinder $r0 $r1 radius $s resolution 30
    graphics $mole sphere $r0 radius $s resolution 30
    graphics $mole sphere $r1 radius $s resolution 30
}

proc triNormal {r0 r1 r2} {
    set a [vecSub $r1 $r0]
    set b [vecSub $r2 $r0]
    set n [vecUnit [vecCross $a $b]]
    return $n
}

proc spherePos {phi theta} {
    set cp [expr {cos($phi)}]
    set sp [expr {sin($phi)}]
    set ct [expr {cos($theta)}]
    set st [expr {sin($theta)}]

    return [list [expr {$ct*$sp}] [expr {$st*$sp}] $cp]
}

proc setColor {color {mole top}} {
    graphics $mole color $color
}

proc clear {{mole top}} {
    graphics $mole delete all
}

proc setMaterial {material color {mole top}} {
    graphics $mole color $color
    graphics $mole material $material
}


proc drawSpheres {pointVar rad {mole top}} {
    upvar $pointVar point

    foreach r $point {
	graphics $mole sphere $r radius $rad resolution 20
    }
}
