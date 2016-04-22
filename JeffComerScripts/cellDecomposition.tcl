# Perform a cell decomposition on a set of points.
# Author: Jeff Comer <jcomer2@illinois.edu>

# Load points in the format "x y z" from a file.
proc loadPoints {fileName} {
    set in [open $fileName r]
    
    set r {}
    foreach line [split [read $in] "\n"] {
	set tok [split [string trim $line]]
	if {[llength $tok] >= 3} {
	    lappend r [list [lindex $tok 0] [lindex $tok 1] [lindex $tok 2]]
	}
    }
    
    close $in
    return $r
}

# Determine the range of the points.
proc getBounds {pos} {
    set xMin [lindex $pos 0 0]
    set xMax [lindex $pos 0 0]
    set yMin [lindex $pos 0 1]
    set yMax [lindex $pos 0 1]
    set zMin [lindex $pos 0 2]
    set zMax [lindex $pos 0 2]

    foreach r $pos {
	foreach {x y z} $r {break}
	
	if {[expr $x < $xMin]} {set xMin $x}
	if {[expr $y < $yMin]} {set yMin $y}
	if {[expr $z < $zMin]} {set zMin $z}
	
	if {[expr $x > $xMax]} {set xMax $x}
	if {[expr $y > $yMax]} {set yMax $y}
	if {[expr $z > $zMax]} {set zMax $z}
    }
    return [list [list $xMin $yMin $zMin] [list $xMax $yMax $zMax]]
}

# Obtain the neighbors' indices for a cell.
proc neighbors {j grid} {
    foreach {lx ly lz} $grid {break}
    set nCells [expr $lx*$ly*$lz]

    set n {}
    for {set kx -1} {$kx <= 1} {incr kx} {
	for {set ky -1} {$ky <= 1} {incr ky} {
	    for {set kz -1} {$kz <= 1} {incr kz} {
		set index [expr $j + $kz + $ky*$lz + $kx*$lz*$ly]
		if {$index >= 0 && $index < $nCells} {
		    lappend n $index
		}
	    }
	}
    }
    
    return $n
}

# Determine the cell index of a point.
proc cellLookup {pos cutoff origin grid} {
    foreach {px py pz} $pos {break}
    foreach {ox oy oz} $origin {break}
    foreach {gx gy gz} $grid {break}

    set ix [expr int(floor(($px-$ox)/$cutoff))]
    set iy [expr int(floor(($py-$oy)/$cutoff))]
    set iz [expr int(floor(($pz-$oz)/$cutoff))]
    return [expr $iz + $iy*$gz + $ix*$gz*$gy]
}

# Obtain the neighbors of a cell.
proc cellNeighbors {grid} {
    set nCells [expr [lindex $grid 0]*[lindex $grid 1]*[lindex $grid 2]]

    # Form the neighbor list.
    set neigh {}
    for {set j 0} {$j < $nCells} {incr j} {
	lappend neigh [neighbors $j $grid] 
    }
    return $neigh
}

# Return the grid origin and grid size.
proc cellGrid {pos cutoff} {
    set bounds [getBounds $pos]
    set dx [expr [lindex $bounds 1 0]-[lindex $bounds 0 0]]
    set dy [expr [lindex $bounds 1 1]-[lindex $bounds 0 1]]
    set dz [expr [lindex $bounds 1 2]-[lindex $bounds 0 2]]
    
    # Determine the origin.
    set ox [expr [lindex $bounds 0 0]-$cutoff]
    set oy [expr [lindex $bounds 0 1]-$cutoff]
    set oz [expr [lindex $bounds 0 2]-$cutoff]

    # Determine the grid size.
    set lx [expr int(ceil($dx/$cutoff))+2]
    set ly [expr int(ceil($dy/$cutoff))+2]
    set lz [expr int(ceil($dz/$cutoff))+2]

    return [list [list $ox $oy $oz] [list $lx $ly $lz]]
}

# Distribute the points into cells.
proc cellDecompose {pos cutoff origin grid cellVar} {
    upvar $cellVar cell
    set nCells [expr [lindex $grid 0]*[lindex $grid 1]*[lindex $grid 2]]

    for {set j 0} {$j < $nCells} {incr j} {
	set cell($j) {}
    }

    foreach r $pos {
	set index [cellLookup $r $cutoff $origin $grid]
	lappend cell($index) $r
    }
}


# Distribute the indices of the points into cells.
proc cellDecomposeIndices {pos cutoff origin grid cellVar} {
    upvar $cellVar cell
    set nCells [expr [lindex $grid 0]*[lindex $grid 1]*[lindex $grid 2]]

    for {set j 0} {$j < $nCells} {incr j} {
	set cell($j) {}
    }

    set i 0
    foreach r $pos {
	set index [cellLookup $r $cutoff $origin $grid]
	lappend cell($index) $i 
	incr i
    }
}

# Count the points in a cell decomposition.
proc cellCountPoints {grid cellVar} {
    upvar $cellVar cell
    set nCells [expr [lindex $grid 0]*[lindex $grid 1]*[lindex $grid 2]]
 
    set n 0
    for {set i 0} {$i < $nCells} {incr i} {
	foreach p $cell($i) {
	    incr n
	}
    }

    return $n
}

# Get the points in the cell decomposition.
proc cellExtractPoints {grid cellVar} {
    upvar $cellVar cell
    set nCells [expr [lindex $grid 0]*[lindex $grid 1]*[lindex $grid 2]]
 
    set points {}
    for {set i 0} {$i < $nCells} {incr i} {
	foreach p $cell($i) {
	    lappend points $p
	}
    }

    return $points
}



