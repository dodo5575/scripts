# Render a solid clip plane.
# to use: vmd -e solidClipPlane.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

proc writeSolidClipPlane {clipid repid molid selText resolution prefix} {
	set center [mol clipplane center $clipid $repid $molid]
	set normal [mol clipplane normal $clipid $repid $molid]
		
	set m [createClipSurf $center $normal $molid $selText $resolution]
	writeMesh m ${prefix}${resolution}.mesh
}

proc solidClipPlane {clipid repid molid selText resolution} {
	set center [mol clipplane center $clipid $repid $molid]
	set normal [mol clipplane normal $clipid $repid $molid]
	set color [mol clipplane color $clipid $repid $molid]
	
	set mesh [createClipSurf $center $normal $molid $selText $resolution]
	drawMesh mesh $molid $color
	return $mesh
}

proc writeClipSurf {center normal molid selText resolution prefix} {
	set m [createClipSurf $center $normal $molid $selText $resolution]
	writeMesh m ${prefix}${resolution}.mesh
}

proc createClipSurf {center normal molid selText resolution} {
	# Parameters
	set distance 4.0
	set offset 0.2
	
	# Prepare the geometry.	
	set normal [vecnorm $normal]	
	set center [vecsub $center [vecscale $offset $normal]]
	set up {0.0 0.0 1.0}
	set v [vecsub $up [vecscale [vecdot $up $normal] $normal]]
		
	# Use another definition of `v' if `normal' is along `up'.
	if {[expr [veclength2 $v] == 0.0]} {
		set up {1.0 0.0 0.0}
		set v [vecsub $up [vecscale [vecdot $up $normal] $normal]]
	}
	puts "Creating a solid clip plane with center $center and normal $normal"
	
	# Define the plane basis.
	set v [vecnorm $v]
	set u [veccross $v $normal]
	puts "The plane basis vector u: $u"
	puts "The plane basis vector v: $v"
	
	# Make more restricted selection text for finding bounds.
	foreach {cx cy cz} $center {nx ny nz} $normal {break}
	set selTextRestr "(${selText}) and (x-$cx)*$nx + (y-$cy)*$ny  + (z-$cz)*$nz < $distance"
	
	# Find the extension of the selection along the plane basis.
	set sel [atomselect $molid $selTextRestr]
	set rAtom [$sel get {x y z}]
	$sel delete
	set ub [getBound rAtom $u $center]
	set vb [getBound rAtom $v $center]
	foreach {u0 u1} $ub {v0 v1} $vb {break}
	puts "Selection extent along u: $u0 to $u1"
	puts "Selection extent along v: $v0 to $v1"
	
	# Make `resolution' triangles along the shorter direction.
	set origin [vecadd $center [vecscale $u0 $u] [vecscale $v0 $v]]
	set du [expr $u1 - $u0]
	set dv [expr $v1 - $v0]
	if {[expr $du < $dv]} {
		set a [expr $du/$resolution]
	} else {
		set a [expr $dv/$resolution]
	}
	set nu [expr int($du/$a) + 1]
	set nv [expr int($dv/($a*0.5*sqrt(3.0))) + 1]
	
	set mesh [createMesh $a $origin $nu $u $nv $v]
	puts "Number of triangles created: [llength $mesh]"
	set t0 [clock seconds]
	set mesh [shadowMesh mesh $molid $selText $distance]
	puts "Shadow mesh computed in [expr [clock seconds] - $t0] s"
	puts "Number of triangles after cutting: [llength $mesh]"
		
	return $mesh
}

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



