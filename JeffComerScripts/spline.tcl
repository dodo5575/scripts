# This script provides procedures for mapping a chain molecule
# along the z-axis to a spline.
# Author: Jeff Comer <jcomer2@illinois.edu>


# Generate a random number on a sphere distributed within angle "angle" of vector "axis".
# This function will take a long time if the angle is very small.
proc sphericalRandom {axis angle} {
	if {$angle <= 0.0} {return $axis0}
	
	set comp [expr cos($angle)]
	# Normalize the axis vector.
	set axis [vecnorm $axis]
	
	# Create a random vector within a half-unit cube and project it onto a sphere.
	set ret [vecnorm [list [expr rand()-0.5] [expr rand()-0.5] [expr rand()-0.5]]]
	
	# Keep the vector within the specified angle of axis0.
	while {[vecdot $ret $axis] < $comp} {
		set ret [vecnorm [list [expr rand()-0.5] [expr rand()-0.5] [expr rand()-0.5]]]
	}
	return $ret
}

# Return a random set of points in a box of side length "side".
proc randomSet {points side} {
	set ret {}
	
	for {set i 0} {$i < $points} {incr i} {
		set r [list [expr rand()-0.5] [expr rand()-0.5] [expr rand()-0.5]]
		set r [vecscale $side $r]
		lappend ret $r
	}
	return $ret
}

# Create a random path of vectors.
proc randomPath {nodes angle} {
	set dr {0 0 1}
	set vert [list {0 0 0}]
	
	for {set i 1} {$i < $nodes} {incr i} {
		# Set a new course "dr" within "angle" of the old "dr".
		set dr [sphericalRandom $dr $angle]
		set r [vecadd [lindex $vert [expr $i-1]] $dr]
		
		# Add the new vertex to the vector list.
		lappend vert $r
		
	}
	return $vert
}

# Return a path forming the letters vmd.
proc VMDPath {} {
return [list\
[list {0 0 0} {-2 5 2} {5 10 0} {8 5 2} {10 0 0}]\
[list {20 10 0} {21 0 0} {25 5 2} {29 0 0} {30 10 0}]\
[list {40 0 0} {40 10 0} {50 7 0} {50 3 0} {44 0 2} {42 -3 5}]]
}

# Add two control points between each vertex.
# Organizes the returned list into sequences of four points with internal vertices repeated.
# Example: r0 c0a c0b r1; r1 c1a c1b r2; r2 c2a c2b r3...
proc insertControlPoints { vert } {
	set last [expr [llength $vert]-1]
	set ret {}
	
	# Set the control point spacing parameter.
	set alpha [expr 1.0/3.0]
	
	# Determine the control points neighboring the internal vertices.		
	for {set i 1} {$i < $last} {incr i} {
		# Determine the control points.
		set ra [lindex $vert [expr $i-1]]
		set r [lindex $vert $i]
		set rb [lindex $vert [expr $i+1]]
		
		set amag [veclength [vecsub $ra $r]]
		set bmag [veclength [vecsub $rb $r]]		
		set aunit [vecnorm [vecsub $ra $r]]
		set bunit [vecnorm [vecsub $rb $r]]
		set side [vecnorm [vecsub $bunit $aunit]]
		
		# These are the control points that neighbor vertex[i].
		set ctrl_a [vecadd $r [vecscale $side [expr -$alpha*$amag]]]
		set ctrl_b [vecadd $r [vecscale $side [expr $alpha*$bmag]]]
		
		# Add the new points to "ret".
		lappend ret $ctrl_a $r $r $ctrl_b
	}
	
	# Set the first control point linearly.
	set r [lindex $vert 0]
	set dmag [veclength [vecsub [lindex $ret 0] $r]] 
	set side [vecnorm [vecsub [lindex $ret 0] $r]]
	set ctrl_b [vecadd $r [vecscale $side [expr $alpha*$dmag]]]
	set ret [concat [list $r] [list $ctrl_b] $ret]
	
	# Set the last control point linearly.
	set r [lindex $vert $last]
	set d [vecsub $r [lindex $ret [expr [llength $ret]-1]]]
	set dmag [veclength $d]
	set side [vecnorm $d]
	set ctrl_a [vecadd $r [vecscale $side [expr -$alpha*$dmag]]]
	lappend ret $ctrl_a $r
	
	# Return the vertex/control point list.
	return $ret
}

# Compute the Bezier spline vectors for each vertex.
# Return sequences of four vectors a, b, c and p0. 
proc bezierParameters { point } {
	set ret {}
	
	foreach {p0 p1 p2 p3} $point {
		set c [vecscale [vecsub $p1 $p0] 3.0]
		set b [vecsub [vecscale [vecsub $p2 $p1] 3.0] $c]
		set a [vecadd $p3 [vecinvert $p0] [vecinvert $c] [vecinvert $b]]
		
		lappend ret $a $b $c $p0
	}
	return $ret
}

# Generate points on a spline given Bezier parameters of vertices.
proc sampleSpline { samples param } {
	set ret {}
	
	foreach {a b c p0} $param {
		for {set i 0} {$i < $samples} {incr i} {
			set t [expr 1.0*$i/$samples]
			set p1 [vecscale $c $t]
			set p2 [vecscale $b [expr $t*$t]]
			set p3 [vecscale $a [expr $t*$t*$t]]
			
			# Locate the point along the spline.
			set r [vecadd $p0 $p1 $p2 $p3]
			
			lappend ret $r
		}
	}
	return $ret	
}

# Create a list recording distance along a series of points.
proc pathLength { curve } {
	set ret {}
	set r0 [lindex $curve 0]
	set s 0
		
	foreach r $curve {
		set s [expr [veclength [vecsub $r $r0]] + $s]
		lappend ret $s
		set r0 $r
	}
	return $ret
}

# Return the total length of the path.
proc totalPathLength { curve } {
	set r0 [lindex $curve 0]
	set s 0
		
	foreach r $curve {
		set s [expr [veclength [vecsub $r $r0]] + $s]
		set r0 $r
	}
	return $s
}

# Create an basis0/axis/angle (matrix/vector/float) list for a series of points.
proc generateBases { curve } {
	set last [expr [llength $curve]-2]
	set ret {}
	
	# Initialize the rotation basis.
	set a [vecnorm [vecsub [lindex $curve 0] [lindex $curve 1]]]
	set basis [transvec $a]
	set basis [transmult $basis [transaxis y 0.5 pi]]
		
	for {set i 0} {$i < $last} {incr i} {
		set r0 [lindex $curve $i]
		set r1 [lindex $curve [expr $i+1]]
		set r2 [lindex $curve [expr $i+2]]
		
		set d0 [vecsub $r1 $r0]
		set d1 [vecsub $r2 $r1]
		
		# Determine the how the basis rotates between the points.
		set comp [vecdot [vecnorm $d0] [vecnorm $d1]]
		if {$comp >= 1.0} {
			set theta 0.0
			set axis [list 1.0 0.0 0.0]
		} else {
			set theta [expr acos($comp)]
			set axis [vecnorm [veccross $d0 $d1]]
		}
		
		lappend ret $basis $axis $theta
				
		# Rotate to the new basis.
		set basis [transmult [transabout $axis $theta rad] $basis]
	}
	
	# Add the tail rotation.
	lappend ret $basis $axis $theta $basis $axis $theta
	return $ret	
}

# Map a list points to a curve given by a list of points.
proc mapToCurve { point curve } {
	set size [llength $point]
	set samples [llength $curve]
	
	# Get the bases at each sample on the curve.
	set arclength [pathLength $curve]
	set rotation [generateBases $curve]
	
	# Find the minimum z-value and the position of the center in the xy-plane.
	set ox 0.0
	set oy 0.0
	set oz [lindex $point 0 2]
	foreach v $point {
		foreach {x y z} $v {break}
		set ox [expr $ox + $x]
		set oy [expr $oy + $y]
		if {$z < $oz} {set oz $z}
	}
	set ox [expr $ox/$size]
	set oy [expr $oy/$size]
	puts [format "curve origin: %g %g %g" $ox $oy $oz]
	
	# Map each point to a curve.
	set p 0
	foreach v $point {
		foreach {x y z} $v {break}
		
		# Calculate coordinates relative to the curve origin.
		set sx [expr $x-$ox]
		set sy [expr $y-$oy]
		set sz [expr $z-$oz]
		
		# Find (sequential search) the samples between which the point lies.
		for {set i 1} {$i < $samples} {incr i} { 
			set s [lindex $arclength $i]
			if {$s > $sz} break
		}
		
		# Check if the point is off the end of the path.
		if {$i == $samples} {continue}
		
		# The previous sample is the basis.
		incr i -1
		set s0 [lindex $arclength $i]
		
		# Calculate the relative progress from the previous sample towards the next.
		set sz [expr $sz - $s0]
		if {[expr $s - $s0] == 0.0} { set t 0.0
		} else {set t [expr $sz/($s - $s0)]}
		
		# Obtain the necessary transformation parameters.
		set basis0 [lindex $rotation [expr 3*$i]]
		set axis [lindex $rotation [expr 3*$i+1]]
		set theta [lindex $rotation [expr 3*$i+2]]
		set pos0 [lindex $curve $i]
		set pos1 [lindex $curve [expr $i+1]]
				
		# Rotate the required part of theta from the initial basis.
		set rot [transabout $axis [expr $theta*$t] rad]
		set rot [transmult $rot $basis0]
		
		# Rotate the point.
		set pos [vectrans $rot [list $sx $sy 0]]
		#set pos [veczero]
				
		# Displace it to the appropriate place on the curve.
		set r [vecadd [vecscale $pos0 [expr 1.0-$t]] [vecscale $pos1 $t]]
		set pos [vecadd $pos $r]
		
		# Modify the point list.
		lset point $p $pos
		incr p	
	}
	return $point
}

# Show the spline by creating a chain of atoms.
proc traceSpline { molid path } {
	set param [bezierParameters [insertControlPoints $path]]
	set spline [sampleSpline 10 $param]
	
	# Check that a sufficient number of atoms are available.
	set all [atomselect $molid "all"]
	if {[$all num] < [llength $spline]} {
		puts "traceSpline failed: Too few atoms"
		return
	}
	
	set last [expr [llength $spline]-1]
	set piece [atomselect $molid [format "index 0 to %d" $last]]
	$piece set {x y z} $spline
}

# Return the basis vectors of a transformation matrix.
proc getBasisX { matrix } {
	return [vectrans $matrix {1.0 0.0 0.0}]
}
proc getBasisY { matrix } {
	return [vectrans $matrix {0.0 1.0 0.0}]
}
proc getBasisZ { matrix } {
	return [vectrans $matrix {0.0 0.0 1.0}]
}

# Show the bases by placing atoms along the axes.
proc traceBases { molid path size } {
	set param [bezierParameters [insertControlPoints $path]]
	set spline [sampleSpline 10 $param]
	set rotation [generateBases $spline]
	set samples [llength $spline]
	
	# Check that a sufficient number of atoms are available.
	set all [atomselect $molid "all"]
	if {[$all num] < [expr 5*[llength $spline]]} {
		puts "traceBases failed: Too few atoms"
		return
	}
	
	set point {}
	for {set i 0} {$i < $samples} {incr i} {
		set basis [lindex $rotation [expr 3*$i]]
		set pos [lindex $spline $i]
		lappend point $pos
		lappend point [vecadd $pos [vecscale $size [getBasisX $basis]]]
		lappend point [vecadd $pos [vecscale $size [getBasisY $basis]]]
		lappend point [vecadd $pos [vecscale -$size [getBasisX $basis]]]
		lappend point [vecadd $pos [vecscale -$size [getBasisY $basis]]]
	}
		
	set last [expr [llength $point]-1]
	set piece [atomselect $molid [format "index 0 to %d" $last]]
	$piece set {x y z} $point
}

proc scalePath { path factor } {
	set ret {}

	foreach r $path {
		lappend ret [vecscale $r $factor]
	}
	return $ret

}

# Scale a path so that its arclength is near that of the z-range of "point".
proc scalePathToZRange { path point } {
	# Find the minimum and maximum z-value of the atoms.
	set zmin [lindex $point 0 2]
	set zmax [lindex $point 0 2]
	foreach v $point {
		set z [lindex $v 2]
		if {$z < $zmin} {set zmin $z}
		if {$z > $zmax} {set zmax $z}
	}
	
	# Rescale with a 1% buffer.
	set arclength [totalPathLength $path]
	set factor [expr ($zmax - $zmin)/$arclength]
	return [scalePath $path [expr $factor*1.01]]
}

proc centerMol { selection } {
	# Get the masses and coordinates of the atoms.
	set data [$selection get {mass x y z}]
	
	# Find the center of mass.
	set mass 0.0
	set center [veczero]
	foreach el $data {
		set m [lindex $el 0]
		set r [list [lindex $el 1] [lindex $el 2] [lindex $el 3]]
		set mass [expr $mass + $m]
		set center [vecadd $center [vecscale $m $r]]
	}
	set center [vecscale [expr 1.0/$mass] $center]
	puts [format "moving center of mass %s to the origin" [list $center]]
	
	# Place the center of mass at the origin.
	$selection moveby [vecinvert $center]
}

# Map a selection of atoms to a spline.
proc mapMolToSpline { selection path } {
	# Get the coordinates of the atoms.
	set point [$selection get {x y z}]
	
	# Generate a spline through the path.
	set param [bezierParameters [insertControlPoints $path]]
	set spline [sampleSpline 40 $param]
		
	# Rescale the spline to match the length of the molecule.
	set spline [scalePathToZRange $spline $point] 
	
	# Perform the mapping.
	set point [mapToCurve $point $spline]
	$selection set {x y z} $point
}




