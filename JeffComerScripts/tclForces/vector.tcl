## namespace exports vector procs used by forces.tcl
namespace eval ::tclForces::vector {
    proc veclength { v } {
	set length2 0
	foreach x $v { lappend length2 "$x*$x" }
	set length2 [join $length2 "+"]
	set length2 [expr $length2]
	return [expr {sqrt($length2)}]
    }
    proc vecnorm { v } {
	#print "veclength $v = [veclength $v]"
	set vl [veclength $v]
	if { $vl == 0 } {print "tclForces: WARNING: vecnorm taken of null vector; returning $v"; return $v}
	#print "sqrt\(vecdot {$v} {$v}\) = [expr sqrt([vecdot $v $v])]"
	return [vecscale [expr {1./$vl}] $v]
    }
    proc vecdot {v1 v2} {
	foreach x1 $v1 x2 $v2 {
	    lappend product "$x1*$x2"
	}
	return [expr [join $product "+"]]
    }
    proc veccross {v1 v2} {
	foreach {x1 x2 x3} $v1 { break }
	foreach {y1 y2 y3} $v2 { break }
	return "[expr {$x2*$y3-$x3*$y2}] [expr {-$x1*$y3+$x3*$y1}] [expr {$x1*$y2-$x2*$y1}]"
    }
    proc veccross2D { v } {
	foreach {x y} $v { break }
	return "[expr {-$y}] $x"
    }
    proc vecprojectSAFE { v axis } {
	set axis [vecnorm $axis]
	return [vecscale [vecdot $v $axis] $axis]
    }
    proc vecproject { v axis } {
	## UNSAFE because it assumes axis is a unit vector
	return [vecscale [vecdot $v $axis] $axis]
    }
    proc vecorthogonal { v axis } {
	return [vecsub $v [vecproject $v $axis]]
    }
    ## List procs
    proc lsum { list } {
	set ret 0
	foreach x $list { lappend ret $x }
	return [expr [join $ret "+"]]
    }

    ## also some math
    proc sign {x} {
	return [expr {($x > 0) - ($x < 0)}]
    }
    
    namespace export *
}


