# Draw lots of things.
# Author: Jeff Comer <jcomer2@illinois.edu>

proc drawFile {fileName} {
    set pos [readPoints $fileName]
    drawSpheres top pos 0.2
}

proc readPoints {fileName} {
    set in [open $fileName r]
    set ret {}

    foreach line [split [read $in] \n] {
	if {[string equal [string index $line 0] "\#"]} {continue}
	set tok [concat $line]
	if {[llength $tok] < 3} {continue}
	lappend ret [lrange $tok 0 2]
    }
    return $ret
}

proc readForces {fileName} {
    set in [open $fileName r]
    set pos {}
    set force {}

    foreach line [split [read $in] \n] {
	if {[string equal [string index $line 0] "\#"]} {continue}
	set tok [concat $line]
	if {[llength $tok] < 6} {continue}
	lappend pos [lrange $tok 0 2]
	lappend force [lrange $tok 3 5]
    }
    return [list $pos $force]
}

proc readField {fileName pointVar forceVar} {
    upvar $pointVar point
    upvar $forceVar force

    set both [readForces $fileName]
    set point [lindex $both 0]
    set force [lindex $both 1]
    return
}

proc drawClear {mole} {
    graphics $mole delete all
}

proc drawSpheres {mole pointVar rad} {
    upvar $pointVar point
    
    foreach r $point {
	graphics $mole sphere $r radius $rad resolution 20
    }
}

proc drawPoints {mole pointVar} {
    upvar $pointVar point
    
    foreach r $point {
	graphics $mole point $r
    }
}

proc drawArrows {mole pointVar forceVar scale} {
    upvar $pointVar point
    upvar $forceVar force

    foreach r $point f $force {
	if {[llength $scale] == 1} {set scale [list $scale $scale $scale]}
	set fMag [veclength $f]
	if {[expr $fMag <= 0.0]} {continue}

	set v [vecscale [expr 1.0/$fMag] $f]
	set v1 [vecscale [expr pow($fMag,0.5)] $scale]
	set arr [list [expr [lindex $v1 0]*[lindex $v 0]] [expr [lindex $v1 1]*[lindex $v 1]] [expr [lindex $v1 2]*[lindex $v 2]]]
	set len [veclength $arr]

	set head [vecadd $r [vecscale 0.6 $arr]]
	set tip [vecadd $head $arr]
	graphics $mole cylinder $head $r radius [expr 0.1*abs($len)] resolution 10
	graphics $mole cone $head $tip radius [expr 0.3*abs($len)] resolution 10
    }
}

proc drawCartoonArrows {mole pointVar forceVar scale view {push 0}} {
    upvar $pointVar point
    upvar $forceVar force

    foreach r $point f $force {
	if {[llength $scale] == 1} {set scale [list $scale $scale $scale]}
	set fMag [veclength $f]
	if {$fMag <= 0.0} {continue}

	set v [vecscale [expr 1.0/$fMag] $f]
	
	# A hack!
	#if {abs([vecdot $v {0 1 0}]) > 0.05} {continue}
	#if {abs([vecdot $v {0 0 1}]) > 0.98} {continue}

	set v1 [vecscale $fMag $scale]
	#set v1 [vecscale [expr pow($fMag,0.5)] $scale]
	set arr [list [expr [lindex $v1 0]*[lindex $v 0]] [expr [lindex $v1 1]*[lindex $v 1]] [expr [lindex $v1 2]*[lindex $v 2]]]
	set len [veclength $arr]

	set head [vecadd $r [vecscale 0.6 $arr]]
	set tip [vecadd $head $arr]
	if {$push == 1} {
	    # Put the tip at the force center.
	    set d [vecsub $r $tip]
	    set r [vecadd $r $d]
	    set head [vecadd $head $d]
	    set tip [vecadd $tip $d]
	} elseif {$push == 2} {
	    # Put the head at the force center.
	    set d [vecsub $r $head]
	    set r [vecadd $r $d]
	    set head [vecadd $head $d]
	    set tip [vecadd $tip $d]
	}

	set radial [veccross $arr $view]
	if {[veclength $radial] <= 0.0} {continue}
	set radial [vecscale [expr 1.0/[veclength $radial]] $radial]
	set radialBase [vecscale [expr 0.15*abs($len)] $radial]
	set radialTip [vecscale [expr 0.4*abs($len)] $radial]
	
	set base0 [vecadd $r $radialBase]
	set base1 [vecadd $r [vecinvert $radialBase]]
	set head0 [vecadd $head $radialBase]
	set head1 [vecadd $head [vecinvert $radialBase]]
	set flank0 [vecadd $head $radialTip]
	set flank1 [vecadd $head [vecinvert $radialTip]]

	# Draw triangles.
	#graphics $mole material Diffuse
	graphics $mole color white
	#graphics $mole triangle $base0 $base1 $head1
	#graphics $mole triangle $head1 $head0 $base0
	#graphics $mole triangle $flank0 $flank1 $tip
	graphics $mole cylinder $head $r radius [expr 0.15*abs($len)] resolution 20
	graphics $mole cone $head $tip radius [expr 0.4*abs($len)] resolution 20
	# Draw lines.
	#set mole 4
	graphics $mole color black
	drawRod $mole $base0 $base1
	drawRod $mole $base1 $head1
	drawRod $mole $head0 $base0
	drawRod $mole $flank0 $head0
	drawRod $mole $head1 $flank1
	drawRod $mole $flank1 $tip
	drawRod $mole $tip $flank0
	#set mole 6
	
    }
}

proc drawRod {mole r0 r1} {
    #graphics $mole line $r0 $r1 width 4
    #graphics $mole material Opaque
    graphics $mole cylinder $r0 $r1 radius 0.3
}

proc setMaterial {mole material color} {
    graphics $mole color $color
    graphics $mole material $material
}

proc setColor {mole color} {
    graphics $mole color $color
}

proc drawGridLines {mole dl p0 p1} {
    set d [vecsub $p1 $p0]
    set nx [expr int(floor([lindex $d 0]/$dl))]
    set ny [expr int(floor([lindex $d 1]/$dl))]
    set nz [expr int(floor([lindex $d 2]/$dl))]

    set rad [expr $dl*0.05]
    set rad1 [expr $dl*0.1]

    for {set ix 0} {$ix < $nx} {incr ix} {
	for {set iy 0} {$iy < $ny} {incr iy} {
	    for {set iz 0} {$iz < $nz} {incr iz} {
		set r0 [vecadd $p0 [list [expr $ix*$dl] [expr $iy*$dl] [expr $iz*$dl]]]
		set r1 [vecadd $r0 [list $dl 0 0]]
		set r2 [vecadd $r0 [list 0 $dl 0]]
		set r3 [vecadd $r0 [list 0 0 $dl]]

		graphics $mole line $r0 $r1 width 2
		graphics $mole line $r0 $r2 width 2
		graphics $mole line $r0 $r3 width 2
	    }
	}
    }
}

proc drawGrid {mole dl p0 p1} {
    set d [vecsub $p1 $p0]
    set nx [expr int(floor([lindex $d 0]/$dl))]
    set ny [expr int(floor([lindex $d 1]/$dl))]
    set nz [expr int(floor([lindex $d 2]/$dl))]

    set rad [expr $dl*0.05]
    set rad1 [expr $dl*0.1]
    
    for {set ix 0} {$ix < $nx} {incr ix} {
	for {set iy 0} {$iy < $ny} {incr iy} {
	    for {set iz 0} {$iz < $nz} {incr iz} {
		set r0 [vecadd $p0 [list [expr $ix*$dl] [expr $iy*$dl] [expr $iz*$dl]]]
		set r1 [vecadd $r0 [list $dl 0 0]]
		set r2 [vecadd $r0 [list 0 $dl 0]]
		set r3 [vecadd $r0 [list 0 0 $dl]]

		#graphics $mole line $r0 $r1 width 2
		#graphics $mole line $r0 $r2 width 2
		#graphics $mole line $r0 $r3 width 2
		
		graphics $mole sphere $r0 radius $rad1 resolution 20
		graphics $mole cylinder $r0 $r1 radius $rad resolution 10
		graphics $mole cylinder $r0 $r2 radius $rad resolution 10
		graphics $mole cylinder $r0 $r3 radius $rad resolution 10
	    }
	}
    }
}

proc drawCells {mole dl p0 p1} {
    set d [vecsub $p1 $p0]
    set nx [expr int(floor([lindex $d 0]/$dl))]
    set ny [expr int(floor([lindex $d 1]/$dl))]
    set nz [expr int(floor([lindex $d 2]/$dl))]

    set rad [expr $dl*0.05]
    set rad1 [expr $dl*0.1]

    set cubeList {{0 1} {0 2} {0 3} {3 5} {3 6} {1 4} {1 5} {5 7} {2 6} {2 4} {4 7} {6 7}}

    for {set ix 0} {$ix < $nx} {incr ix} {
	for {set iy 0} {$iy < $ny} {incr iy} {
	    for {set iz 0} {$iz < $nz} {incr iz} {
		set r0 [vecadd $p0 [list [expr $ix*$dl] [expr $iy*$dl] [expr $iz*$dl]]]
		set r(0) $r0
		set r(1) [vecadd $r0 [list $dl 0 0]]
		set r(2) [vecadd $r0 [list 0 $dl 0]]
		set r(3) [vecadd $r0 [list 0 0 $dl]]
		set r(4) [vecadd $r0 [list $dl $dl 0]]
		set r(5) [vecadd $r0 [list $dl 0 $dl]]
		set r(6) [vecadd $r0 [list 0 $dl $dl]]
		set r(7) [vecadd $r0 [list $dl $dl $dl]]

		for {set i 0} {$i < 8} {incr i} {
		    graphics $mole sphere $r($i) radius $rad1 resolution 20
		}
		
		foreach pair $cubeList {
		    foreach {vert0 vert1} $pair {break}
		    graphics $mole cylinder $r($vert0) $r($vert1) radius $rad resolution 10
		}
	    }
	}
    }
}



