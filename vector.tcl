# Vector procedures
# All procedures work on n-vectors except
# vecCross, matInvert, vecZero, matIdentity, matMake4, matRandomRot
# which assume 3-vectors.
# Author: Jeff Comer <jcomer2@illinois.edu>

proc vecZero {} {
    return [list 0.0 0.0 0.0]
}
proc vecInvert {a} {
    set b {}
    foreach ai $a {
	lappend b [expr -$ai]
    }
    return $b
}
proc vecAdd {a b} {
    set c {}
    foreach ai $a bi $b {
	lappend c [expr $ai+$bi]
    }
    return $c
}
proc vecSub {a b} {
    set c {}
    foreach ai $a bi $b {
	lappend c [expr $ai-$bi]
    }
    return $c
}
proc vecDot {a b} {
    set sum 0
    foreach ai $a bi $b {
	set sum [expr $sum + $ai*$bi]
    }
    return $sum
}
proc vecCross {a b} {
    foreach {ax ay az} $a {break}
    foreach {bx by bz} $b {break}
    
    set cx [expr $ay*$bz - $az*$by]
    set cy [expr $az*$bx - $ax*$bz]
    set cz [expr $ax*$by - $ay*$bx]
    return [list $cx $cy $cz]
}
proc vecLength {a} {
    set sum 0
    foreach ai $a {
	set sum [expr $sum + $ai*$ai]
    }
    return [expr sqrt($sum)]
}
proc vecLength2 {a} {
    set sum 0
    foreach ai $a {
	set sum [expr $sum + $ai*$ai]
    }
    return [expr $sum]
}
proc vecScale {s a} {
    set b {}
    foreach ai $a {
	lappend b [expr $ai*$s]
    }
    return $b
}
proc vecTransform {m a} {
    set b {}
    foreach row $m {
	lappend b [vecDot $row $a]
    }
    return $b
}
proc matIdentity {} {
    return [list [list 1.0 0.0 0.0] [list 0.0 1.0 0.0] [list 0.0 0.0 1.0]]
}
proc matTranspose {m} {
    set n [llength $m]
    set t $m

    for {set i 0} {$i < $n} {incr i} {
	for {set j 0} {$j < $n} {incr j} {
	    lset t $i $j [lindex $m $j $i]
	}
    }
    return $t
}
proc matScale {s m} {
    set ret {}
    
    foreach row $m {
	set newRow {}
	foreach n $row {
	    lappend newRow [expr $s*$n]
	}
	lappend ret $newRow
    }
    return $ret
}

proc matDet {m} {

    foreach {mxx mxy mxz} [lindex $m 0] {break}
    foreach {myx myy myz} [lindex $m 1] {break}
    foreach {mzx mzy mzz} [lindex $m 2] {break}
    
    set det [expr 1.0*($mxx*($myy*$mzz-$myz*$mzy) - $mxy*($myx*$mzz-$myz*$mzx) + $mxz*($myx*$mzy-$myy*$mzx))]

    return $det
}

proc matInvert {m} {
    foreach {mxx mxy mxz} [lindex $m 0] {break}
    foreach {myx myy myz} [lindex $m 1] {break}
    foreach {mzx mzy mzz} [lindex $m 2] {break}
    
    set det [expr 1.0*($mxx*($myy*$mzz-$myz*$mzy) - $mxy*($myx*$mzz-$myz*$mzx) + $mxz*($myx*$mzy-$myy*$mzx))]
    set ixx [expr ($myy*$mzz - $myz*$mzy)/$det]
    set ixy [expr -($mxy*$mzz - $mxz*$mzy)/$det]
    set ixz [expr ($mxy*$myz - $mxz*$myy)/$det]
    set iyx [expr -($myx*$mzz - $myz*$mzx)/$det]
    set iyy [expr ($mxx*$mzz - $mxz*$mzx)/$det]
    set iyz [expr -($mxx*$myz - $mxz*$myx)/$det]
    set izx [expr ($myx*$mzy - $myy*$mzx)/$det]
    set izy [expr -($mxx*$mzy - $mxy*$mzx)/$det]
    set izz [expr ($mxx*$myy - $mxy*$myx)/$det]

    return [list [list $ixx $ixy $ixz] [list $iyx $iyy $iyz] [list $izx $izy $izz]]
}
proc matMul {a b} {
    set bt [matTranspose $b]

    set ret {}
    foreach rowA $a {
	set r {}
	foreach colB $bt {
	    lappend r [vecDot $rowA $colB]
	}
	lappend ret $r
    }
    return $ret
}
proc matRandomRot {} {
    set pi [expr 4.*atan(1.)]
    # Create a random rotation matrix, uniform on a sphere.
    set a [expr 2.0*rand()*$pi]
    set b [expr acos(2.0*rand()-1.0)]
    set c [expr 2.0*rand()*$pi]

    set ca [expr cos($a)]
    set sa [expr sin($a)]
    set za [expr -sin($a)]
    set cb [expr cos($b)]
    set sb [expr sin($b)]
    set zb [expr -sin($b)]
    set cc [expr cos($c)]
    set sc [expr sin($c)]
    set zc [expr -sin($c)]

    set ta [list [list $ca $za 0] [list $sa $ca 0] [list 0 0 1]]
    set tb [list [list 1 0 0] [list 0 $cb $zb] [list 0 $sb $cb]]
    set tc [list [list $cc $zc 0] [list $sc $cc 0] [list 0 0 1]]

    set basis $ta
    set basis [matMul $tb $basis]
    set basis [matMul $tc $basis]
    return $basis
}
proc matMake4 {m {d {0.0 0.0 0.0}}} {
    set ret {}
    lappend ret [concat [lindex $m 0] [lindex $d 0]]
    lappend ret [concat [lindex $m 1] [lindex $d 1]]
    lappend ret [concat [lindex $m 2] [lindex $d 2]]
    lappend ret [list 0.0 0.0 0.0 1.0]
}



