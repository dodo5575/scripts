# Render a solid clip plane.
# to use: vmd -e solidClipPlane.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

source vector.tcl

# Read a mesh from a file.
proc readMesh {fileName} {
    set in [open $fileName r]
    set mesh {}	
    foreach line [split [read $in] \n] {
	set item [split [string trim $line]]
	if {[llength $item] == 9} {
	    set p0 [list [lindex $item 0] [lindex $item 1] [lindex $item 2]]
	    set p1 [list [lindex $item 3] [lindex $item 4] [lindex $item 5]]
	    set p2 [list [lindex $item 6] [lindex $item 7] [lindex $item 8]]
	    lappend mesh [list $p0 $p1 $p2]
	}
    }
    close $in
    return $mesh
}

# Read a mesh from a file.
proc readQuads {fileName} {
    set in [open $fileName r]
    set mesh {}	
    foreach line [split [read $in] \n] {
	set item [split [string trim $line]]
	if {[llength $item] == 12} {
	    set p0 [lrange $item 0 2]
	    set p1 [lrange $item 3 5]
	    set p2 [lrange $item 6 8]
	    set p3 [lrange $item 9 11]
	    set quad [list $p0 $p1 $p2 $p3]

	    foreach tri [decomposeQuad $quad] {
		lappend mesh $tri
	    }
	}
	if {[llength $item] == 9} {
	    set p0 [lrange $item 0 2]
	    set p1 [lrange $item 3 5]
	    set p2 [lrange $item 6 8]
	    lappend mesh [list $p0 $p1 $p2]
	}
    }
    close $in
    return $mesh
}

# Decompose a quad into four triangles.
proc decomposeQuad {quad} {
    set ret {}

    set c [veczero]
    foreach p $quad {
	set c [vecadd $c $p]
    }
    set c [vecscale 0.25 $c]

    lappend ret [list [lindex $quad 0] [lindex $quad 1] $c]
    lappend ret [list [lindex $quad 1] [lindex $quad 2] $c]
    lappend ret [list [lindex $quad 2] [lindex $quad 3] $c]
    lappend ret [list [lindex $quad 3] [lindex $quad 0] $c]
    return $ret
}

# Write a mesh to a file.
proc writeMesh {meshVar fileName} {
    upvar $meshVar mesh
    
    set out [open $fileName w]
    foreach tri $mesh {
	foreach {p0 p1 p2} $tri {break}
	puts $out "$p0 $p1 $p2"
    }
    close $out
}

# Find the bound of the points along a vector.
proc getBound {pointVar v center} {
    upvar $pointVar point
    set v0 [vecdot $v [vecsub [lindex $point 0] $center]]
    set v1 $v0
    
    foreach r $point {
	set vi [vecdot $v [vecsub $r $center]]

	if {[expr $vi < $v0]} {set v0 $vi}
	if {[expr $vi > $v1]} {set v1 $vi}
    }
    
    return [list $v0 $v1]
}

# Create a mesh of triangles in a plane.
proc createMesh {a origin nu u nv v} {
    set vScale [expr 0.5*sqrt(3.0)]
    set col [vecscale $a $u]
    set row [vecadd [vecscale [expr $vScale*$a] $v] [vecscale 0.5 $col]] 
    
    set mesh {}
    for {set j 0} {$j < $nv} {incr j} {
	# Shift every other row.
	if {[expr $j%2]} {
	    set shift -0.5
	} else {
	    set shift 0.0
	}
	
	# Determine the positions of the triangle vertices.
	for {set i 0} {$i < $nu} {incr i} {
	    set x [vecscale [expr $a*($i+$shift)] $u]
	    set y [vecscale [expr $a*$j*$vScale] $v]
	    set p0 [vecadd $x $y $origin]
	    set p1 [vecadd $p0 $col]
	    set p2 [vecadd $p0 $row]
	    set p3 [vecadd $p2 $col]
	    lappend mesh [list $p0 $p1 $p2]
	    lappend mesh [list $p1 $p3 $p2]
	}	
    }
    
    return $mesh
}


# Cut out triangles that are not within `distance' of a selected atom.
proc shadowMesh {meshVar molid selText distance} {
    upvar $meshVar mesh
    set third [expr 1.0/3.0]
    
    set shadow {}
    foreach tri $mesh {
	foreach {p0 p1 p2} $tri {break}
	set center [vecscale $third [vecadd $p0 $p1 $p2]]
	
	if {[nearAtoms $center $molid $selText $distance] > 0} {
	    lappend shadow $tri
	}
	
    }
    
    return $shadow
}


# How many atoms are within `distance' of point?
proc nearAtoms {pos molid selText distance} {
    foreach {x y z} $pos {break}
    set dSq [expr $distance*$distance]
    set sel [atomselect $molid "($selText) and (x-$x)^2 + (y-$y)^2 + (z-$z)^2 < $dSq"]
    set n [$sel num]
    $sel delete
    
    return $n
}

# Draw a mesh using VMD's `graphics' commands.
proc drawVertices {meshVar molid color} {
    upvar $meshVar mesh
    
    graphics $molid color $color
    foreach tri $mesh {
	foreach {p0 p1 p2} $tri {break}
	graphics $molid sphere $p0
	graphics $molid sphere $p1
	graphics $molid sphere $p2
    }
}

# Draw a mesh using VMD's `graphics' commands.
proc drawMesh {meshVar molid color} {
    upvar $meshVar mesh
    
    graphics $molid color $color
    foreach tri $mesh {
	foreach {p0 p1 p2} $tri {break}
	graphics $molid triangle $p0 $p1 $p2
    }
}

# Draw a mesh using VMD's `graphics' commands.
proc drawMeshTransform {meshVar molid color origin basis} {
    upvar $meshVar mesh
    
    graphics $molid color $color
    foreach tri $mesh {
	foreach {p0 p1 p2} $tri {break}
	set r0 [vecAdd [vecTransform $basis $p0] $origin]
	set r1 [vecAdd [vecTransform $basis $p1] $origin]
	set r2 [vecAdd [vecTransform $basis $p2] $origin]

	graphics $molid triangle $r0 $r1 $r2
    }
}
