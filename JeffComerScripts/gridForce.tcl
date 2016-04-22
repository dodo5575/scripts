# gridForce.tcl
# Grid force procedures
# Many procedures use a grid structure with the variables:
#  nx, ny, nz # nodes along each grid axis
#  origin # grid origin 3-vector in real space
#  delta # grid basis 3x3 matrix
#  data # list of potential values z fast, y medium, and x slow.
#
#  # Extra variables for convenience
#  size # nx*ny*nz
#  deltaInv # inverse grid basis 3x3 matrix
#  
# Author: Jeff Comer <jcomer2@illinois.edu>

# Returns the potential value at each point in gridVar(data),
# the number of grid points along each
# cell axis in gridVar(nx,ny,nz), the lattice vectors in
# gridVar(delta) (3x3 matrix), and the cell origin in
# gridVar(origin) (3-vector).
# Data order is z fast, y medium, and x slow.

proc newGridDim {gridVar gridX gridY gridZ lx ly lz} {
    upvar $gridVar grid

    set grid(data) {}
    set grid(origin) [vecScale -0.5 [list $lx $ly $lz]]
    set grid(nx) [expr int(ceil($lx/$gridX))]
    set grid(ny) [expr int(ceil($ly/$gridY))]
    set grid(nz) [expr int(ceil($lz/$gridZ))]
    set grid(size) [expr $grid(nx)*$grid(ny)*$grid(nz)]
    set grid(delta) [list [list $gridX 0 0] [list 0 $gridY 0] [list 0 0 $gridZ]]
    set grid(deltaInv) [matInvert $grid(delta)]

    for {set i 0} {$i < $grid(size)} {incr i} {
	lappend grid(data) 0
    }
    return $grid(size)
}

proc newGrid {gridVar delta nx ny nz} {
    upvar $gridVar grid
    
    set grid(data) {}
    set grid(origin) [vecScale -0.5 [vecTransform $delta [list $nx $ny $nz]]]
    set grid(nx) $nx
    set grid(ny) $ny
    set grid(nz) $nz
    set grid(size) [expr $grid(nx)*$grid(ny)*$grid(nz)]
    set grid(delta) $delta
    set grid(deltaInv) [matInvert $grid(delta)]

    for {set i 0} {$i < $grid(size)} {incr i} {
	lappend grid(data) 0
    }
    return $grid(size)
}

proc newGridFit {gridVar lx ly lz dl} {
    upvar $gridVar grid

    set nx [expr int(ceil($lx/$dl))]
    set ny [expr int(ceil($ly/$dl))]
    set nz [expr int(ceil($lz/$dl))]
    set gridX [expr double($lx)/$nx]
    set gridY [expr double($ly)/$ny]
    set gridZ [expr double($lz)/$nz]
    
    set grid(data) {}
    set grid(origin) [vecScale -0.5 [list $lx $ly $lz]]
    set grid(nx) $nx
    set grid(ny) $ny
    set grid(nz) $nz
    set grid(size) [expr $grid(nx)*$grid(ny)*$grid(nz)]
    set grid(delta) [list [list $gridX 0 0] [list 0 $gridY 0] [list 0 0 $gridZ]]
    set grid(deltaInv) [matInvert $grid(delta)]

    for {set i 0} {$i < $grid(size)} {incr i} {
	lappend grid(data) 0.0
    }
    return $grid(size)
}

proc newGridFitXYZ {gridVar lx ly lz dx dy dz} {
    upvar $gridVar grid

    set nx [expr int(ceil($lx/$dx))]
    set ny [expr int(ceil($ly/$dy))]
    set nz [expr int(ceil($lz/$dz))]
    set gridX [expr double($lx)/$nx]
    set gridY [expr double($ly)/$ny]
    set gridZ [expr double($lz)/$nz]
    
    set grid(data) {}
    set grid(origin) [vecScale -0.5 [list $lx $ly $lz]]
    set grid(nx) $nx
    set grid(ny) $ny
    set grid(nz) $nz
    set grid(size) [expr $grid(nx)*$grid(ny)*$grid(nz)]
    set grid(delta) [list [list $gridX 0 0] [list 0 $gridY 0] [list 0 0 $gridZ]]
    set grid(deltaInv) [matInvert $grid(delta)]

    for {set i 0} {$i < $grid(size)} {incr i} {
	lappend grid(data) 0.0
    }
    return $grid(size)
}

proc newGridBox {gridVar ax ay az dl} {
    upvar $gridVar grid

    set dx [vecLength $ax]
    set dy [vecLength $ay]
    set dz [vecLength $az]

    set nx [expr int(ceil($dx/$dl))]
    set ny [expr int(ceil($dy/$dl))]
    set nz [expr int(ceil($dz/$dl))]
    
    set ex [vecScale [expr 1.0/$nx] $ax]
    set ey [vecScale [expr 1.0/$ny] $ay]
    set ez [vecScale [expr 1.0/$nz] $az]
   
    set grid(data) {}
    set grid(nx) $nx
    set grid(ny) $ny
    set grid(nz) $nz
    set grid(size) [expr $grid(nx)*$grid(ny)*$grid(nz)]
    set grid(delta) [matTranspose [list $ex $ey $ez]]
    set grid(deltaInv) [matInvert $grid(delta)]
    set diag [vecTransform $grid(delta) [list $nx $ny $nz]]
    set grid(origin) [vecScale -0.5 $diag]

    for {set i 0} {$i < $grid(size)} {incr i} {
	lappend grid(data) 0.0
    }
    return $grid(size)
}

proc newGridXsc {gridVar xscFile dl} {
    upvar $gridVar grid

    # Read the system size from the xsc file.
    set in [open $xscFile r]
    foreach line [split [read $in] "\n"] {
	if {![string match "#*" $line]} {
	    set param [split $line]
	    puts $param
	    set a [lrange $param 1 3]
	    set b [lrange $param 4 6]
	    set c [lrange $param 7 9]
	    break
	}
    }
    close $in

    newGridBox grid $a $b $c $dl
}

proc newGridCuboid {gridVar r0 r1 dx dy dz} {
    upvar $gridVar grid
    
    set l [vecSub $r1 $r0]
    foreach {lx ly lz} $l { break }
    set nx [expr int(ceil($lx/$dx))]
    set ny [expr int(ceil($ly/$dy))]
    set nz [expr int(ceil($lz/$dz))]
    set gridX [expr double($lx)/$nx]
    set gridY [expr double($ly)/$ny]
    set gridZ [expr double($lz)/$nz]
    
    set grid(data) {}
    set grid(origin) $r0
    set grid(nx) $nx
    set grid(ny) $ny
    set grid(nz) $nz
    set grid(size) [expr $grid(nx)*$grid(ny)*$grid(nz)]
    set grid(delta) [list [list $gridX 0 0] [list 0 $gridY 0] [list 0 0 $gridZ]]
    set grid(deltaInv) [matInvert $grid(delta)]

    for {set i 0} {$i < $grid(size)} {incr i} {
	lappend grid(data) 0.0
    }
    return $grid(size)
}

proc copyGridDim {gridVar newVar} {
    upvar $newVar new
    upvar $gridVar grid

    set new(origin) $grid(origin)
    set new(nx) $grid(nx)
    set new(ny) $grid(ny)
    set new(nz) $grid(nz)
    set new(size) $grid(size)
    set new(delta) $grid(delta)
    set new(deltaInv) $grid(deltaInv)
    zeroGrid new
    return $new(size)
}

proc zeroGrid {gridVar} {
    upvar $gridVar grid
    
    set grid(data) {}
    for {set i 0} {$i < $grid(size)} {incr i} {
	lappend grid(data) 0.0
    }
}

proc readDx {gridVar fileName} {
    upvar $gridVar grid
    set in [open $fileName r]

    set items -1
    set n -1
    set grid(data) {}
    set grid(delta) {}

    foreach line [split [read $in] "\n"] {
	set tok [concat $line]
	
	# Ignore comments and blank lines.
	if {[string match "\#*" $tok]} {continue}
	if {[string length $line] < 1} {continue}
	if {$items > 0 && $items == $n} {continue}

	if {$n < $items} {
	    # Read a line of numeric data.
	    foreach t $tok {
		lappend grid(data) $t
		incr n
	    }
	} elseif {[string match "object*" $line]} {
	    # Read an object.
	    if {[llength $tok] < 4} {
		puts "Warning: Invalid dx object."
	    }
	    
	    set type [lindex $tok 3]
	    if {[string equal $type "gridpositions"]} {
		# Read the number of nodes along the grid axes.
		set grid(nx) [lindex $tok 5]
		set grid(ny) [lindex $tok 6]
		set grid(nz) [lindex $tok 7]
	    } elseif {[string equal $type "array"]} {
		# Begin reading the numeric data.
		set items [lindex $tok 9]
		set n 0
		set grid(size) $items
		if {[expr $grid(nx)*$grid(ny)*$grid(nz) != $grid(size)]} {
		    puts "Warning: Number of items $grid(size) not consistent with grid dimensions $grid(nx) $grid(ny) $grid(nz)"
		}
	    }
	} elseif {[string match "origin*" $line]} {
	    # Read the grid origin.
	    set grid(origin) [lrange $tok 1 3]
	    
	} elseif {[string match "delta*" $line]} {
	    # Read the grid basis.
	    lappend grid(delta) [lrange $tok 1 3]
	}

    }
    close $in

    # Transpose the delta matrix.
    set grid(delta) [matTranspose $grid(delta)]
    # Set the inverse delta matrix.
    set grid(deltaInv) [matInvert $grid(delta)]

    return $grid(size)
}

proc readDxHeader {gridVar fileName} {
    upvar $gridVar grid
    set in [open $fileName r]

    set items -1
    set n -1
    set grid(data) {}
    set grid(delta) {}

    foreach line [split [read $in] "\n"] {
	set tok [concat $line]
	
	# Ignore comments and blank lines.
	if {[string match "\#*" $tok]} {continue}
	if {[string length $line] < 1} {continue}
	if {$items > 0 && $items == $n} {continue}

	if {$n < $items} {
	    # Quit if we reach a data line.
	    break
	} elseif {[string match "object*" $line]} {
	    # Read an object.
	    if {[llength $tok] < 4} {
		puts "Warning: Invalid dx object."
	    }
	    
	    set type [lindex $tok 3]
	    if {[string equal $type "gridpositions"]} {
		# Read the number of nodes along the grid axes.
		set grid(nx) [lindex $tok 5]
		set grid(ny) [lindex $tok 6]
		set grid(nz) [lindex $tok 7]
	    } elseif {[string equal $type "array"]} {
		# Begin reading the numeric data.
		set items [lindex $tok 9]
		set n 0
		set grid(size) $items
		if {[expr $grid(nx)*$grid(ny)*$grid(nz) != $grid(size)]} {
		    puts "Warning: Number of items $grid(size) not consistent with grid dimensions $grid(nx) $grid(ny) $grid(nz)"
		}
	    }
	} elseif {[string match "origin*" $line]} {
	    # Read the grid origin.
	    set grid(origin) [lrange $tok 1 3]
	    
	} elseif {[string match "delta*" $line]} {
	    # Read the grid basis.
	    lappend grid(delta) [lrange $tok 1 3]
	}

    }
    close $in

    # Transpose the delta matrix.
    set grid(delta) [matTranspose $grid(delta)]
    # Set the inverse delta matrix.
    set grid(deltaInv) [matInvert $grid(delta)]

    return $grid(size)
}


# Writes a dx file given a grid structure as returned by readDx.
proc writeDx {gridVar fileName} {
    upvar $gridVar grid

    foreach {ox oy oz} $grid(origin) {break}
    set delta [join $grid(delta)]
    foreach {exx eyx ezx exy eyy ezy exz eyz ezz} $delta {break}

    set out [open $fileName w]
    # Write the headers.
    puts $out "\# NAMD gridforce grid"
    # Write the grid attributes.
    puts $out "object 1 class gridpositions counts $grid(nx) $grid(ny) $grid(nz)"
    puts $out "origin $ox $oy $oz"
    puts $out "delta $exx $exy $exz"
    puts $out "delta $eyx $eyy $eyz"
    puts $out "delta $ezx $ezy $ezz"
    puts $out "object 2 class gridconnections counts $grid(nx) $grid(ny) $grid(nz)"
    puts $out "object 3 class array type double rank 0 items $grid(size) data follows"
    
    # Write the data.
    set penultima [expr 3*($grid(size)/3)]
    set mod [expr $grid(size)-$penultima]
    set i 0
    foreach {v0 v1 v2} $grid(data) {
	puts $out "$v0 $v1 $v2"
	
	incr i 3
	# Treat the last line specially.
	if {$i >= $penultima} {break}
    }

    if {$mod == 1} {
	puts $out "[lindex $grid(data) end]"
    } elseif {$mod == 2} {
	puts $out "[lindex $grid(data) end-1] [lindex $grid(data) end]"
    }

    close $out
}

# Writes the dx file header given an output stream.
proc writeDxHeader {gridVar out} {
    upvar $gridVar grid

    foreach {ox oy oz} $grid(origin) {break}
    set delta [join $grid(delta)]
    foreach {exx eyx ezx exy eyy ezy exz eyz ezz} $delta {break}

    # Write the headers.
    puts $out "\# NAMD gridforce grid"
    # Write the grid attributes.
    puts $out "object 1 class gridpositions counts $grid(nx) $grid(ny) $grid(nz)"
    puts $out "origin $ox $oy $oz"
    puts $out "delta $exx $exy $exz"
    puts $out "delta $eyx $eyy $eyz"
    puts $out "delta $ezx $ezy $ezz"
    puts $out "object 2 class gridconnections counts $grid(nx) $grid(ny) $grid(nz)"
    puts $out "object 3 class array type double rank 0 items $grid(size) data follows"
}

# Get the extent of the grid in world space.
proc getGridSize {gridVar} {
    upvar $gridVar grid
    set r [vecTransform $grid(delta) [list $grid(nx) $grid(ny) $grid(nz)]]
    return $r
}

# Get the length of a lattice cell diagonal in world space.
proc getVoxelSize {gridVar} {
    upvar $gridVar grid
    set r [vecTransform $grid(delta) {1 1 1}]
    return $r
}

proc getVolume {gridVar} {
    upvar $gridVar grid
    return [expr {abs([matDet $grid(delta)])}]
}

# Get the basis matrix for the entire system.
proc getSystemCell {gridVar} {
    upvar $gridVar grid
    set s [list [list $grid(nx) 0.0 0.0] [list 0.0 $grid(ny) 0.0] [list 0.0 0.0 $grid(nz)]]
    return [matMul $s $grid(delta)]
}

# Convert from lattice space to world space.
proc getPosition {gridVar rl} {
    upvar $gridVar grid
    set r [vecAdd [vecTransform $grid(delta) $rl] $grid(origin)]
    return $r
}

# Convert from world space to lattice space.
proc getLatticePosition {gridVar r} {
    upvar $gridVar grid
    set ri [vecTransform $grid(deltaInv) [vecSub $r $grid(origin)]]
    return $ri
}

# Get the lattice space position of the home (using floor) node.
proc homeLatticePosition {gridVar pos} {
    upvar $gridVar grid
    set l [vecTransform $grid(deltaInv) [vecSub $pos $grid(origin)]]
    foreach {lx ly lz} $l {break}
 
    # Find the home node.
    set ix [expr {int(floor($lx))}]
    set iy [expr {int(floor($ly))}]
    set iz [expr {int(floor($lz))}]
    return [list $ix $iy $iz]
}

# Get the node index of the home (using floor) node.
proc homeIndex {gridVar pos} {
    upvar $gridVar grid
    set l [vecTransform $grid(deltaInv) [vecSub $pos $grid(origin)]]
    foreach {lx ly lz} $l {break}

    # Find the home node.
    set ix [expr {int(floor($lx))}]
    set iy [expr {int(floor($ly))}]
    set iz [expr {int(floor($lz))}]

    # Wrap on the boundaries.
    set ix [expr {$ix % $grid(nx)}]
    set iy [expr {$iy % $grid(ny)}]
    set iz [expr {$iz % $grid(nz)}]

    return [expr {$iz + $grid(nz)*$iy + $grid(nz)*$grid(ny)*$ix}]
}

# Is the point in the grid?
proc inGrid {gridVar pos} {
    upvar $gridVar grid
    set l [vecTransform $grid(deltaInv) [vecSub $pos $grid(origin)]]
    foreach {lx ly lz} $l {break}

    # Find the home node.
    set ix [expr {int(floor($lx))}]
    set iy [expr {int(floor($ly))}]
    set iz [expr {int(floor($lz))}]

    if {$ix < 0 || $ix >= $grid(nx)} { return 0 }
    if {$iy < 0 || $iy >= $grid(ny)} { return 0 }
    if {$iz < 0 || $iz >= $grid(nz)} { return 0 }
    return 1
}


# Get the node index of the nearest node.
proc nearestIndex {gridVar pos} {
    upvar $gridVar grid
    set l [vecTransform $grid(deltaInv) [vecSub $pos $grid(origin)]]
    foreach {lx ly lz} $l {break}

    # Find the home node.
    set ix [expr {int(floor($lx+0.5))}]
    set iy [expr {int(floor($ly+0.5))}]
    set iz [expr {int(floor($lz+0.5))}]

    # Wrap on the boundaries.
    set ix [expr {$ix % $grid(nx)}]
    set iy [expr {$iy % $grid(ny)}]
    set iz [expr {$iz % $grid(nz)}]

    return [expr {$iz + $grid(nz)*$iy + $grid(nz)*$grid(ny)*$ix}]
}

# Get the lattice space position of the closest node.
proc nearestLatticePosition {gridVar pos} {
    upvar $gridVar grid
    set l [vecTransform $grid(deltaInv) [vecSub $pos $grid(origin)]]
    foreach {lx ly lz} $l {break}
 
    # Find the home node.
    set ix [expr {int(floor($lx+0.5))}]
    set iy [expr {int(floor($ly+0.5))}]
    set iz [expr {int(floor($lz+0.5))}]
    return [list $ix $iy $iz]
}


# Convert a lattice position to a lattice index.
# Returns -1 for positions off the grid.
proc latticeIndex1 {gridVar rl} {
    upvar $gridVar grid
    foreach {ix iy iz} $rl { break }
 
    if {$ix < 0 || $ix >= $grid(nx)} {return -1}
    if {$iy < 0 || $iy >= $grid(ny)} {return -1}
    if {$iz < 0 || $iz >= $grid(nz)} {return -1}
    return [expr {$iz + $grid(nz)*$iy + $grid(nz)*$grid(ny)*$ix}]
}

# Convert a lattice position to a lattice index.
# Returns -1 for positions off the grid.
proc latticeIndex {gridVar ix iy iz} {
    upvar $gridVar grid
    
    if {$ix < 0 || $ix >= $grid(nx)} {return -1}
    if {$iy < 0 || $iy >= $grid(ny)} {return -1}
    if {$iz < 0 || $iz >= $grid(nz)} {return -1}
    return [expr {$iz + $grid(nz)*$iy + $grid(nz)*$grid(ny)*$ix}]
}

# Convert an index to a lattice position.
# Returns -1 for positions off the grid.
proc latticePosition {gridVar index} {
    upvar $gridVar grid

    if {$index < 0 || $index >= $grid(size)} {return -1}
    set iz [expr {$index%$grid(nz)}]
    set iy [expr {($index/$grid(nz))%$grid(ny)}]
    set ix [expr {($index/$grid(nz)/$grid(ny))%$grid(nx)}]
    return [list $ix $iy $iz]
}

proc indexToWorld {gridVar index} {
    upvar $gridVar grid

    if {$index < 0 || $index >= $grid(size)} {return -1}
    set iz [expr {$index%$grid(nz)}]
    set iy [expr {($index/$grid(nz))%$grid(ny)}]
    set ix [expr {($index/$grid(nz)/$grid(ny))%$grid(nx)}]
    set r [vecAdd [vecTransform $grid(delta) [list $ix $iy $iz]] $grid(origin)]
    return $r
}

proc worldToIndex {gridVar pos} {
    upvar $gridVar grid
    set l [vecTransform $grid(deltaInv) [vecSub $pos $grid(origin)]]
    foreach {lx ly lz} $l {break}

    # Find the home node.
    set ix [expr {int(floor($lx+0.5))}]
    set iy [expr {int(floor($ly+0.5))}]
    set iz [expr {int(floor($lz+0.5))}]

    # Wrap on the boundaries.
    set ix [expr {$ix % $grid(nx)}]
    set iy [expr {$iy % $grid(ny)}]
    set iz [expr {$iz % $grid(nz)}]

    return [expr {$iz + $grid(nz)*$iy + $grid(nz)*$grid(ny)*$ix}]
}

proc wrapReal {x l} {
    set l [expr {double($l)}]
    set image [expr {int(floor($x/$l))}]
    set x [expr {$x - $image*$l}]

    return $x
}

proc wrapDiffReal {x l} {
    set l [expr {double($l)}]
    set image [expr {int(floor($x/$l))}]
    set x [expr {$x - $image*$l}]

    if {$x >= 0.5*$l} { set x [expr {$x - $l}] }
    return $x
}

proc wrapDiff {gridVar r} {
    upvar $gridVar grid

    set rl [vecTransform $grid(deltaInv) $r]
    foreach {x y z} $rl { break }
    
    set x [wrapDiffReal $x $grid(nx)]
    set y [wrapDiffReal $y $grid(ny)]
    set z [wrapDiffReal $z $grid(nz)]

    return [vecTransform $grid(delta) [list $x $y $z]]
}

proc wrap {gridVar r} {
    upvar $gridVar grid
    set rl [vecTransform $grid(deltaInv) [vecSub $r $grid(origin)]] 
    foreach {x y z} $rl { break }
    
    set x [wrapReal $x $grid(nx)]
    set y [wrapReal $y $grid(ny)]
    set z [wrapReal $z $grid(nz)]

    return [vecAdd [vecTransform $grid(delta) [list $x $y $z]] $grid(origin)]
}

# Find the grid value of the node closest to a point.
# No interpolation is done.
proc getPotential {gridVar r} {
    upvar $gridVar grid

    set ri [vecTransform $grid(deltaInv) [vecSub $r $grid(origin)]]
    
    set ix [expr {int(floor([lindex $ri 0]+0.5))}]
    set iy [expr {int(floor([lindex $ri 1]+0.5))}]
    set iz [expr {int(floor([lindex $ri 2]+0.5))}]

    set i [expr {$iz + $grid(nz)*$iy + $grid(nz)*$grid(ny)*$ix}]
    if {$i < 0 || $i >= $grid(size)} {return 0.0}
    return [lindex $grid(data) $i]
}

# Sample 'n' points along a line segment from 'r0' to 'r1' to generate
# a one-dimensional profile of the grid.
proc potentialProfile {gridVar r0 r1 n fileName} {
    upvar $gridVar grid
    set dr [vecScale [expr 1./$n] [vecSub $r1 $r0]]
    set d [vecLength $dr]
    set d0 [expr 0.5*$n*$d]
    
    set out [open $fileName w]
    set ret {}
    for {set i 0} {$i < $n} {incr i} {
	set r [vecAdd $r0 [vecScale $i $dr]]
	set pot [interpolatePotential grid $r]
	
	puts $out "[expr $i*$d-$d0] $pot"
    }
    close $out
    return
}

# Sample 'n' points along a line segment from 'r0' to 'r1' to generate
# a one-dimensional profile of the grid.
proc potentialProfileZ {gridVar r0 r1 n fileName} {
    upvar $gridVar grid
    set dr [vecScale [expr 1./$n] [vecSub $r1 $r0]]
    set d [vecLength $dr]
    set d0 [expr 0.5*$n*$d]
    
    set out [open $fileName w]
    set ret {}
    for {set i 0} {$i < $n} {incr i} {
	set r [vecAdd $r0 [vecScale $i $dr]]
	set pot [interpolatePotential grid $r]
	
	puts $out "[lindex $r 2] $pot"
    }
    close $out
    return
}

# Sample 'n' points along a line segment from 'r0' to 'r1' to generate
# a one-dimensional profile of the negative gradient of the grid.
proc forceProfile {gridVar r0 r1 n dir fileName} {
    upvar $gridVar grid
    set dr [vecScale [expr 1./$n] [vecSub $r1 $r0]]
    set d [vecLength $dr]
    set d0 [expr 0.5*$n*$d]
    
    set out [open $fileName w]
    set ret {}
    for {set i 0} {$i < $n} {incr i} {
	set r [vecAdd $r0 [vecScale $i $dr]]
	set pot [interpolateForce grid $r $dir]
	
	puts $out "[expr $i*$d-$d0] $pot"
    }
    close $out
    return
}

# Sample 'n*n' points in a plane to generate
# a two-dimensional profile of the grid.
proc potentialCrossSection {gridVar lx ly} {
}

# Decrease the resolution of the grid by an integer factor 'step'.
# No interpolation is done.
proc coarsenGrid {gridVar newVar step} {
    upvar $gridVar grid
    upvar $newVar new

    set step [expr int(floor($step))]
    set new(nx) [expr $grid(nx)/$step]
    set new(ny) [expr $grid(ny)/$step]
    set new(nz) [expr $grid(nz)/$step]
    set new(size) [expr $new(nx)*$new(ny)*$new(nz)]
    set new(origin) $grid(origin)
    set new(delta) [matScale $step $grid(delta)]
    set new(deltaInv) [matInvert $new(delta)]

    set new(data) {}
    for {set ix 0} {$ix < $new(nx)} {incr ix} {
	for {set iy 0} {$iy < $new(ny)} {incr iy} {
	    for {set iz 0} {$iz < $new(nz)} {incr iz} {
		set j [expr {$step*($iz + $iy*$grid(nz) + $ix*$grid(ny)*$grid(nz))}]
		lappend new(data) [lindex $grid(data) $j]
	    }
	}
    }
    return
}

# Add two grids, interpolating the second onto the first.
# The new grid will have the dimensions of the first.
proc addGrids {gridVar addVar newVar} {
    upvar $gridVar grid
    upvar $addVar add
    upvar $newVar new

    copyGridDim grid new

    set j 0
    set new(data) {}
    for {set ix 0} {$ix < $grid(nx)} {incr ix} {
	for {set iy 0} {$iy < $grid(ny)} {incr iy} {
	    for {set iz 0} {$iz < $grid(nz)} {incr iz} {
		set r [getPosition new [list $ix $iy $iz]]
		#set pot1 [interpolatePotential add $r]
		set pot1 [interpolatePotentialLinear add $r]
		set pot [expr {[lindex $grid(data) $j] + $pot1}]
		lappend new(data) $pot
		incr j
	    }
	}
    }   
    return
}

# Add two grids with the same dimensions.
proc addGridsSame {gridVar addVar newVar} {
    upvar $gridVar grid
    upvar $addVar add
    upvar $newVar new

    copyGridDim grid new

    set j 0
    set new(data) {}
    for {set ix 0} {$ix < $grid(nx)} {incr ix} {
	for {set iy 0} {$iy < $grid(ny)} {incr iy} {
	    for {set iz 0} {$iz < $grid(nz)} {incr iz} {
		set pot [expr {[lindex $grid(data) $j] + [lindex $add(data) $j]}]
		lappend new(data) $pot
		incr j
	    }
	}
    }
    return
}


# Trim off zeros from grid.
proc trimGrid {gridVar newVar {tol 0.0}} {
    upvar $gridVar grid
    upvar $newVar new

    copyGridDim grid new

    set ix0 -1
    set ix1 0
    set iy0 -1
    set iy1 0
    set iz0 -1
    set iz1 0

    # Find the trimming boundaries.
    set j 0
    for {set ix 0} {$ix < $grid(nx)} {incr ix} {
	for {set iy 0} {$iy < $grid(ny)} {incr iy} {
	    for {set iz 0} {$iz < $grid(nz)} {incr iz} {
		set pot [lindex $grid(data) $j]
		
		if {abs($pot) > $tol} {
		    if {$ix0 < 0} {set ix0 $ix}
		    set ix1 $ix

		    if {$iy0 < 0} {set iy0 $iy}
		    set iy1 $iy
		    
		    if {$iz0 < 0} {set iz0 $iz}
		    set iz1 $iz
		}
		incr j
	    }
	}
    }

    if {$ix0 < 0 || $iy0 < 0 || $iz0 < 0} {
	set new(nx) 0
	set new(ny) 0
	set new(nz) 0
	set new(size) 0
	set new(data) {}
	return
    }

    # Trim the grid.
    set new(data) {}
    for {set ix $ix0} {$ix <= $ix1} {incr ix} {
	for {set iy $iy0} {$iy <= $iy1} {incr iy} {
	    for {set iz $iz0} {$iz <= $iz1} {incr iz} {
		set j [expr {$iz + $iy*$grid(nz) + $ix*$grid(ny)*$grid(nz)}]
      		lappend new(data) [lindex $grid(data) $j]
	    }
	}
    }

    # Set the new dimensions.
    set new(nx) [expr $ix1 - $ix0 + 1]
    set new(ny) [expr $iy1 - $iy0 + 1]
    set new(nz) [expr $iz1 - $iz0 + 1]
    set new(size) [expr $new(nx)*$new(ny)*$new(nz)]
    set dl [list $ix0 $iy0 $iz0]
    set new(origin) [vecAdd $grid(origin) [vecTransform $grid(delta) $dl]]

    #puts "initial size: $grid(nx) $grid(ny) $grid(nz)"
    #puts "final size: $new(nx) $new(ny) $new(nz)"
    return
}


# Resample a grid to an orthogonal lattice with interpolation.
proc samplePotential {gridVar newVar nx ny nz} {
    upvar $gridVar grid
    upvar $newVar new

    set new(nx) $nx
    set new(ny) $ny
    set new(nz) $nz
    set new(size) [expr $nx*$ny*$nz]
    set new(origin) $grid(origin)
    
    set last [list [expr $grid(nx)-1] [expr $grid(ny)-1] [expr $grid(nz)-1]]
    set d1 [vecTransform $grid(delta) $last]
    set dx [expr [lindex $d1 0]/$nx]
    set dy [expr [lindex $d1 1]/$ny]
    set dz [expr [lindex $d1 2]/$nz]

    set new(data) {}
    for {set ix 0} {$ix < $nx} {incr ix} {
	for {set iy 0} {$iy < $ny} {incr iy} {
	    for {set iz 0} {$iz < $nz} {incr iz} {
		set r [list [expr {$dx*$ix}] [expr {$dy*$iy}] [expr {$dz*$iz}]]
		set r [vecAdd $r $new(origin)]
		lappend new(data) [interpolatePotential grid $r]
		#lappend new(data) [getPotential grid $r]
	    }
	}
    }

    set new(delta) [list [list $dx 0 0] [list 0 $dy 0] [list 0 0 $dz]]
    # Set the inverse delta matrix.
    set new(deltaInv) [matInvert $new(delta)]
    return
}

# Resample a grid to an orthogonal lattice over a given domain.
proc samplePotentialVolume {gridVar newVar r0 r1 dl} {
    upvar $gridVar grid
    upvar $newVar new

    set d [vecSub $r1 $r0]
    set nx [expr int(floor([lindex $d 0]/$dl))]
    set ny [expr int(floor([lindex $d 1]/$dl))]
    set nz [expr int(floor([lindex $d 2]/$dl))]
    set new(nx) $nx
    set new(ny) $ny
    set new(nz) $nz
    set new(size) [expr $nx*$ny*$nz]
    set new(origin) $r0
   
    set new(data) {}
    for {set ix 0} {$ix < $nx} {incr ix} {
	for {set iy 0} {$iy < $ny} {incr iy} {
	    for {set iz 0} {$iz < $nz} {incr iz} {
		set r [list [expr {$dl*$ix}] [expr {$dl*$iy}] [expr {$dl*$iz}]]
		set r [vecAdd $r $new(origin)]
		lappend new(data) [interpolatePotential grid $r]
	    }
	}
    }

    set new(delta) [list [list $dl 0 0] [list 0 $dl 0] [list 0 0 $dl]]
    # Set the inverse delta matrix.
    set new(deltaInv) [matInvert $new(delta)]
    return
}

# Resample a grid to an orthogonal lattice with interpolation of the
# negative gradient.
proc sampleForce {gridVar newVar nx ny nz dir} {
    upvar $gridVar grid
    upvar $newVar new
    
    set do [vecTransform $grid(delta) [list 2 2 2]] 
    set new(nx) $nx
    set new(ny) $ny
    set new(nz) $nz
    set new(size) [expr $nx*$ny*$nz]
    set new(origin) [vecAdd $grid(origin) $do]
    
    set last [list [expr $grid(nx)-5] [expr $grid(ny)-5] [expr $grid(nz)-5]]
    set d1 [vecTransform $grid(delta) $last]
    set dx [expr [lindex $d1 0]/$nx]
    set dy [expr [lindex $d1 1]/$ny]
    set dz [expr [lindex $d1 2]/$nz]

    set new(data) {}
    for {set ix 0} {$ix < $nx} {incr ix} {
	for {set iy 0} {$iy < $ny} {incr iy} {
	    for {set iz 0} {$iz < $nz} {incr iz} {
		set r [list [expr {$dx*$ix}] [expr {$dy*$iy}] [expr {$dz*$iz}]]
		set r [vecAdd $r $new(origin)]
		lappend new(data) [interpolateForce grid $r $dir]
	    }
	}
    }

    set new(delta) [list [list $dx 0 0] [list 0 $dy 0] [list 0 0 $dz]]
    # Set the inverse delta matrix.
    set new(deltaInv) [matInvert $new(delta)]
    return
}

# Resample a grid to an orthogonal lattice with interpolation of the
# negative gradient.
proc sampleForceVolume {gridVar newVar r0 r1 dl dir} {
    upvar $gridVar grid
    upvar $newVar new

    set d [vecSub $r1 $r0]
    set nx [expr int(floor([lindex $d 0]/$dl))]
    set ny [expr int(floor([lindex $d 1]/$dl))]
    set nz [expr int(floor([lindex $d 2]/$dl))]
    set new(nx) $nx
    set new(ny) $ny
    set new(nz) $nz
    set new(size) [expr $nx*$ny*$nz]
    set new(origin) $r0
   
    set new(data) {}
    for {set ix 0} {$ix < $nx} {incr ix} {
	for {set iy 0} {$iy < $ny} {incr iy} {
	    for {set iz 0} {$iz < $nz} {incr iz} {
		set r [list [expr {$dl*$ix}] [expr {$dl*$iy}] [expr {$dl*$iz}]]
		set r [vecAdd $r $new(origin)]
		lappend new(data) [interpolateForce grid $r $dir]
	    }
	}
    }

    set new(delta) [list [list $dl 0 0] [list 0 $dl 0] [list 0 0 $dl]]
    # Set the inverse delta matrix.
    set new(deltaInv) [matInvert $new(delta)]
    return
}

# Expand a grid to contain additional points in chosen
# directions with the value of the boundary.
# These are useful as ghost points in interpolation schemes.
# The component of the gradient normal to the boundary is
# then zero.
proc ghostGrid {gridVar newVar {ghost {2 2 2}}} {
    upvar $gridVar grid
    upvar $newVar new
    
    # Handle scalar and vector ghost parameters.
    if {[llength $ghost] < 3} {
	set ghostX $ghost
	set ghostY $ghost
	set ghostZ $ghost
    } else {
	foreach {ghostX ghostY ghostZ} $ghost {break}
    }

    set new(delta) $grid(delta)
    set new(deltaInv) $grid(deltaInv)
    set new(nx) [expr $grid(nx)+2*$ghostX]
    set new(ny) [expr $grid(ny)+2*$ghostY]
    set new(nz) [expr $grid(nz)+2*$ghostZ]
    
    set new(data) {}
    for {set ix 0} {$ix < $new(nx)} {incr ix} {
	for {set iy 0} {$iy < $new(ny)} {incr iy} {
	    for {set iz 0} {$iz < $new(nz)} {incr iz} {
		set dx [expr {($ix-$ghostX)}]
		set dy [expr {($iy-$ghostY)}]
		set dz [expr {($iz-$ghostZ)}]

		# Fix the ghost points to the boundary value.
		if {$dx < 0} {
		    set dx 0
		} elseif {$dx >= $grid(nx)} {
		    set dx [expr {$grid(nx)-1}]
		}
		if {$dy < 0} {
		    set dy 0
		} elseif {$dy >= $grid(ny)} {
		    set dy [expr {$grid(ny)-1}]
		}
		if {$dz < 0} {
		    set dz 0
		} elseif {$dz >= $grid(nz)} {
		    set dz [expr {$grid(nz)-1}]
		}
		set i [expr {$dz + $dy*$grid(nz) + $dx*$grid(nz)*$grid(ny)}]
		
		lappend new(data) [lindex $grid(data) $i]
	    }
	}
    }

    set new(size) [llength $new(data)]
    set do [vecTransform $new(delta) [list $ghostX $ghostY $ghostZ]] 
    set new(origin) [vecSub $grid(origin) $do]
    return
}


# Expand a grid to contain additional zero points in chosen
# directions with the value of the boundary.
# Negative pad values can also be used to crop the grid.
proc padGrid {gridVar newVar {pad {2 2 2}}} {
    upvar $gridVar grid
    upvar $newVar new
    
    # Handle scalar and vector pad parameters.
    if {[llength $pad] < 3} {
	set padX $pad
	set padY $pad
	set padZ $pad
    } else {
	foreach {padX padY padZ} $pad {break}
    }

    set new(delta) $grid(delta)
    set new(deltaInv) $grid(deltaInv)
    set new(nx) [expr $grid(nx)+2*$padX]
    set new(ny) [expr $grid(ny)+2*$padY]
    set new(nz) [expr $grid(nz)+2*$padZ]
    
    set new(data) {}
    for {set ix 0} {$ix < $new(nx)} {incr ix} {
	for {set iy 0} {$iy < $new(ny)} {incr iy} {
	    for {set iz 0} {$iz < $new(nz)} {incr iz} {
		set dx [expr {$ix-$padX}]
		set dy [expr {$iy-$padY}]
		set dz [expr {$iz-$padZ}]

		# Fill the pad points with zero.
		if {$dx < 0 || $dx >= $grid(nx) || $dy < 0 || $dy >= $grid(ny)
		    || $dz < 0 || $dz >= $grid(nz)} {
		    lappend new(data) 0.0
		} else {
		    set i [expr {$dz + $dy*$grid(nz) + $dx*$grid(nz)*$grid(ny)}]
		    
		    lappend new(data) [lindex $grid(data) $i]
		}
	    }
	}
    }

    set new(size) [llength $new(data)]
    set do [vecTransform $new(delta) [list $padX $padY $padZ]] 
    set new(origin) [vecSub $grid(origin) $do]
    return
}

proc distributeLinear {gridVar pos {weight 1.0}} {
    upvar $gridVar grid
    set l [vecTransform $grid(deltaInv) [vecSub $pos $grid(origin)]]
    foreach {lx ly lz} $l {break}
 
    # Find the home node.
    set ix [expr int(floor($lx))]
    set iy [expr int(floor($ly))]
    set iz [expr int(floor($lz))]
    if {$ix < 0 || $ix >= [expr $grid(nx)-1]} {return}
    if {$iy < 0 || $iy >= [expr $grid(ny)-1]} {return}
    if {$iz < 0 || $iz >= [expr $grid(nz)-1]} {return}

    # Find the interpolation coordinates.
    set wx(1) [expr $lx - $ix]
    set wy(1) [expr $ly - $iy]
    set wz(1) [expr $lz - $iz]
    set wx(0) [expr 1.0-$wx(1)]
    set wy(0) [expr 1.0-$wy(1)]
    set wz(0) [expr 1.0-$wz(1)]
    
    # Distribute the weight among the nodes.
    for {set i 0} {$i <= 1} {incr i} {
	for {set j 0} {$j <= 1} {incr j} {
	    for {set k 0} {$k <= 1} {incr k} {
		set neigh [expr $k+$iz + ($j+$iy)*$grid(nz) + ($i+$ix)*$grid(nz)*$grid(ny)]
		set w [expr {$weight*$wx($i)*$wy($j)*$wz($k)}]
		set pot [lindex $grid(data) $neigh]
		lset grid(data) $neigh [expr {$pot + $w}]
	    }
	}
    }
}

proc interpolatePotentialLinear {gridVar pos} {
    upvar $gridVar grid
    set l [vecTransform $grid(deltaInv) [vecSub $pos $grid(origin)]]
    foreach {lx ly lz} $l {break}
 
    # Find the home node.
    set ix [expr {int(floor($lx))}]
    set iy [expr {int(floor($ly))}]
    set iz [expr {int(floor($lz))}]
    if {$ix < 0 || $ix >= [expr {$grid(nx)-1}]} {return 0.0}
    if {$iy < 0 || $iy >= [expr {$grid(ny)-1}]} {return 0.0}
    if {$iz < 0 || $iz >= [expr {$grid(nz)-1}]} {return 0.0}

    # Find the interpolation coordinates.
    set wx [expr {$lx - $ix}]
    set wy [expr {$ly - $iy}]
    set wz [expr {$lz - $iz}]

    # Find the function value at the neighbors.
    for {set i 0} {$i <= 1} {incr i} {
	for {set j 0} {$j <= 1} {incr j} {
	    for {set k 0} {$k <= 1} {incr k} {
		set neigh [expr {$k+$iz + ($j+$iy)*$grid(nz) + ($i+$ix)*$grid(nz)*$grid(ny)}]
		set f(${i}${j}${k}) [lindex $grid(data) $neigh]
	    }
	}
    }
    
    # Mix along lattice vector x.
    set g00 [expr {$f(100)*$wx + (1.0-$wx)*$f(000)}]
    set g10 [expr {$f(110)*$wx + (1.0-$wx)*$f(010)}]
    set g01 [expr {$f(101)*$wx + (1.0-$wx)*$f(001)}]
    set g11 [expr {$f(111)*$wx + (1.0-$wx)*$f(011)}]

    # Mix along lattice vector y.
    set h0 [expr {$g10*$wy + (1.0-$wy)*$g00}]
    set h1 [expr {$g11*$wy + (1.0-$wy)*$g01}]

    # Mix along lattice vector z.
    return [expr {$h1*$wz + (1.0-$wz)*$h0}]
}

# Compute the force at a given position.
# A cubic interpolant is used.
proc interpolatePotential {gridVar pos} {
    upvar $gridVar grid
    set l [vecTransform $grid(deltaInv) [vecSub $pos $grid(origin)]]
    foreach {lx ly lz} $l {break}
 
    # Find the home node.
    set ix [expr {int(floor($lx))}]
    set iy [expr {int(floor($ly))}]
    set iz [expr {int(floor($lz))}]
    if {$ix < 1 || $ix >= [expr {$grid(nx)-2}]} {return 0.0}
    if {$iy < 1 || $iy >= [expr {$grid(ny)-2}]} {return 0.0}
    if {$iz < 1 || $iz >= [expr {$grid(nz)-2}]} {return 0.0}

    # Find the interpolation coordinates.
    set wx [expr {$lx - $ix}]
    set wy [expr {$ly - $iy}]
    set wz [expr {$lz - $iz}]

    # Find the function value at the neighbors.
    for {set i 0} {$i < 4} {incr i} {
	for {set j 0} {$j < 4} {incr j} {
	    for {set k 0} {$k < 4} {incr k} {
		set neigh [expr {$k-1+$iz + ($j-1+$iy)*$grid(nz) + ($i-1+$ix)*$grid(nz)*$grid(ny)}]
		set f(${i}${j}${k}) [lindex $grid(data) $neigh]
	    }
	}
    }

    # Mix along x.
    for {set j 0} {$j < 4} {incr j} {
	for {set k 0} {$k < 4} {incr k} {
	    set a3 [expr {0.5*(-$f(0${j}${k}) + 3*$f(1${j}${k}) - 3*$f(2${j}${k}) + $f(3${j}${k}))}]
	    set a2 [expr {0.5*(2*$f(0${j}${k}) - 5*$f(1${j}${k}) + 4*$f(2${j}${k}) - $f(3${j}${k}))}]
	    set a1 [expr {0.5*(-$f(0${j}${k}) + $f(2${j}${k}))}]
	    set a0 [expr {$f(1${j}${k})}]
	
	    set f(${j}${k}) [expr {$a3*$wx*$wx*$wx + $a2*$wx*$wx + $a1*$wx + $a0}]
	}
    }

    # Mix along y.
    for {set k 0} {$k < 4} {incr k} {
	set a3 [expr {0.5*(-$f(0${k}) + 3*$f(1${k}) - 3*$f(2${k}) + $f(3${k}))}]
	set a2 [expr {0.5*(2*$f(0${k}) - 5*$f(1${k}) + 4*$f(2${k}) - $f(3${k}))}]
	set a1 [expr {0.5*(-$f(0${k}) + $f(2${k}))}]
	set a0 [expr {$f(1${k})}]
	
	set f(${k}) [expr {$a3*$wy*$wy*$wy + $a2*$wy*$wy + $a1*$wy + $a0}]
    }

    # Mix along z.
    set a3 [expr {0.5*(-$f(0) + 3*$f(1) - 3*$f(2) + $f(3))}]
    set a2 [expr {0.5*(2*$f(0) - 5*$f(1) + 4*$f(2) - $f(3))}]
    set a1 [expr {0.5*(-$f(0) + $f(2))}]
    set a0 [expr {$f(1)}]

    return [expr {$a3*$wz*$wz*$wz + $a2*$wz*$wz + $a1*$wz + $a0}]
}

# Compute the force at a given position along a given direction.
# An analytic derivative of a cubic interpolant is taken
# in the desired direction, so that the negative gradient
# along that direction is returned.
# dir = 0, 1, 2 gives the x, y, z lattice axes, respectively.
proc interpolateForce {gridVar pos dir} {
    upvar $gridVar grid
    set l [vecTransform $grid(deltaInv) [vecSub $pos $grid(origin)]]
    foreach {lx ly lz} $l {break}
 
    # Find the home node.
    set ix [expr {int(floor($lx))}]
    set iy [expr {int(floor($ly))}]
    set iz [expr {int(floor($lz))}]
    if {$ix < 1 || $ix >= [expr {$grid(nx)-2}]} {return 0.0}
    if {$iy < 1 || $iy >= [expr {$grid(ny)-2}]} {return 0.0}
    if {$iz < 1 || $iz >= [expr {$grid(nz)-2}]} {return 0.0}

    # Find the interpolation coordinates.
    # The derivative is taken along the direction $dir.
    # Swap the axes so that the derivative will be taken along the
    # direction called x.
    if {$dir == 0} {
	set inx $ix
	set iny $iy
	set inz $iz
	set bx [expr {$grid(nz)*$grid(ny)}]
	set by $grid(nz)
	set bz 1
	set wx [expr {$lx - $ix}]
	set wy [expr {$ly - $iy}]
	set wz [expr {$lz - $iz}]
    } elseif {$dir == 1} {
	set inx $iy
	set iny $ix
	set inz $iy
	set bx $grid(nz)
	set by 1
	set bz [expr {$grid(nz)*$grid(ny)}]
	set wx [expr {$ly - $iy}]
	set wy [expr {$lz - $iz}]
	set wz [expr {$lx - $ix}]
    } elseif {$dir == 2} {
	set inx $iz
	set iny $ix
	set inz $iy
	set bx 1
	set by [expr {$grid(nz)*$grid(ny)}]
	set bz $grid(nz)
	set wx [expr {$lz - $iz}]
	set wy [expr {$lx - $ix}]
	set wz [expr {$ly - $iy}]
    } else {
	puts "Error gridForce::interpolateForce: Direction must be 0, 1, or 2."
	return
    }

    # Find the function value at the neighbors.
    for {set i 0} {$i < 4} {incr i} {
	for {set j 0} {$j < 4} {incr j} {
	    for {set k 0} {$k < 4} {incr k} {
		set neigh [expr {($k-1+$inz)*$bz + ($j-1+$iny)*$by + ($i-1+$inx)*$bx}]
		set f(${i}${j}${k}) [lindex $grid(data) $neigh]
	    }
	}
    }

    # Mix along x, taking the derivative.
    for {set j 0} {$j < 4} {incr j} {
	for {set k 0} {$k < 4} {incr k} {
	    set a3 [expr {0.5*(-$f(0${j}${k}) + 3*$f(1${j}${k}) - 3*$f(2${j}${k}) + $f(3${j}${k}))}]
	    set a2 [expr {0.5*(2*$f(0${j}${k}) - 5*$f(1${j}${k}) + 4*$f(2${j}${k}) - $f(3${j}${k}))}]
	    set a1 [expr {0.5*(-$f(0${j}${k}) + $f(2${j}${k}))}]
	    
	    set f(${j}${k}) [expr {3.0*$a3*$wx*$wx + 2.0*$a2*$wx + $a1}]
	}
    }

    # Mix along y.
    for {set k 0} {$k < 4} {incr k} {
	set a3 [expr {0.5*(-$f(0${k}) + 3*$f(1${k}) - 3*$f(2${k}) + $f(3${k}))}]
	set a2 [expr {0.5*(2*$f(0${k}) - 5*$f(1${k}) + 4*$f(2${k}) - $f(3${k}))}]
	set a1 [expr {0.5*(-$f(0${k}) + $f(2${k}))}]
	set a0 [expr {$f(1${k})}]
	
	set f(${k}) [expr {$a3*$wy*$wy*$wy + $a2*$wy*$wy + $a1*$wy + $a0}]
    }

    # Mix along z.
    set a3 [expr {0.5*(-$f(0) + 3*$f(1) - 3*$f(2) + $f(3))}]
    set a2 [expr {0.5*(2*$f(0) - 5*$f(1) + 4*$f(2) - $f(3))}]
    set a1 [expr {0.5*(-$f(0) + $f(2))}]
    set a0 [expr {$f(1)}]

    return [expr {-($a3*$wz*$wz*$wz + $a2*$wz*$wz + $a1*$wz + $a0)}]
}

# Take the negative gradient of a grid using a second
# order derivative.
proc gridForce {potVar forceVarX forceVarY forceVarZ} {
    upvar $potVar pot
    upvar $forceVarX fx
    upvar $forceVarY fy
    upvar $forceVarZ fz

    # Prepare the resulting fields.
    set do [vecTransform $pot(delta) [list 1 1 1]] 
    set fx(nx) [expr $pot(nx)-2]
    set fx(ny) [expr $pot(ny)-2]
    set fx(nz) [expr $pot(nz)-2]
    set fx(delta) $pot(delta)
    set fx(deltaInv) $pot(deltaInv)
    set fx(origin) [vecAdd $pot(origin) $do]
    set fx(size) [expr $fx(nx)*$fx(ny)*$fx(nz)]

    set fy(nx) $fx(nx)
    set fy(nx) $fx(nx)
    set fy(ny) $fx(ny)
    set fy(nz) $fx(nz)
    set fy(delta) $fx(delta)
    set fy(deltaInv) $fx(deltaInv)
    set fy(origin) $fx(origin)
    set fy(size) $fx(size)

    set fz(nx) $fx(nx)
    set fz(nx) $fx(nx)
    set fz(ny) $fx(ny)
    set fz(nz) $fx(nz)
    set fz(delta) $fx(delta)
    set fz(deltaInv) $fx(deltaInv)
    set fz(origin) $fx(origin)
    set fz(size) $fx(size)

    # Compute the derivative along each lattice axis.
    array set hop [list 0 [expr $pot(nz)*$pot(ny)] 1 $pot(nz) 2 1]
    array set f [list 0 {} 1 {} 2 {}]
    
    for {set ix 0} {$ix < $fx(nx)} {incr ix} {
	for {set iy 0} {$iy < $fx(ny)} {incr iy} {
	    for {set iz 0} {$iz < $fx(nz)} {incr iz} {
		set potMe [expr {$iz+1 + ($iy+1)*$pot(nz) + ($ix+1)*$pot(nz)*$pot(ny)}]
		
		for {set b 0} {$b < 3} {incr b} {
		    set potLeft [lindex $pot(data) [expr {$potMe-$hop($b)}]]
		    set potRight [lindex $pot(data) [expr {$potMe+$hop($b)}]]
		    lappend f($b) [expr {-0.5*($potRight - $potLeft)}]
		}
	    }
	}
    }

    # Transform into world coordinates.
    set deltaInvTrans [matTranspose $pot(deltaInv)]
    foreach f0 $f(0) f1 $f(1) f2 $f(2) {
	set gradL [list $f0 $f1 $f2]
	set gradR [vecTransform $deltaInvTrans $gradL]
	lappend fx(data) [lindex $gradR 0]
	lappend fy(data) [lindex $gradR 1]
	lappend fz(data) [lindex $gradR 2]
    }
    return
}

# Given a vector field, return a magnitude field.
proc gridMagnitude {forceVarX forceVarY forceVarZ magVar} {
    upvar $magVar mag
    upvar $forceVarX fx
    upvar $forceVarY fy
    upvar $forceVarZ fz

    set mag(nx) $fx(nx)
    set mag(ny) $fx(ny)
    set mag(nz) $fx(nz)
    set mag(delta) $fx(delta)
    set mag(deltaInv) $fx(deltaInv)
    set mag(origin) $fx(origin)
    set mag(size) $fx(size)

    set mag(data) {}
    foreach x $fx(data) y $fy(data) z $fz(data) {
	set f [expr {sqrt($x*$x + $y*$y + $z*$z)}]
	lappend mag(data) $f
    }
    return
}

proc averageCrossSection {gridVar iz} {
    upvar $gridVar grid

    # Compute the average potential in this slice.
    set pot 0.0
    for {set iy 0} {$iy < $grid(ny)} {incr iy} {
  	for {set ix 0} {$ix < $grid(nx)} {incr ix} {
	    set i [expr {$iz + $grid(nz)*$iy + $grid(nz)*$grid(ny)*$ix}]
	    set pot [expr {$pot + [lindex $grid(data) $i]}]
	}
    }

    set pot [expr $pot/($grid(nx)*$grid(ny))]
    return $pot
}

# Add a constant field to the grid.
proc addConstantPotential {gridVar pot0} {
    upvar $gridVar grid
    
    set newData {}
    foreach pot $grid(data) {
	lappend newData [expr {$pot + $pot0}]
    }

    set grid(data) $newData
    return
}

# Add a constant field to the grid.
proc addConstantPotential {gridVar pot0} {
    upvar $gridVar grid
    
    set newData {}
    foreach pot $grid(data) {
	lappend newData [expr {$pot + $pot0}]
    }

    set grid(data) $newData
    return
}

# Add a linearly varying potential in the direction of eField0.
proc addConstantForce {gridVar eField0} {
    upvar $gridVar grid
    
    if {[llength $eField0] != 3} {
	puts "ERROR:gridForce:addLinearPotential: eField0 must be a 3-vector."
        return
    }
    
    # Add -E0.r to the grid so that the gradient is E0.
    set j 0
    set newData {}
    for {set ix 0} {$ix < $grid(nx)} {incr ix} {
	for {set iy 0} {$iy < $grid(ny)} {incr iy} {
	    for {set iz 0} {$iz < $grid(nz)} {incr iz} {
		set dr [vecTransform $grid(delta) [list $ix $iy $iz]]
		set pot0 [expr -[vecDot $dr $eField0]]
		set pot [lindex $grid(data) $j]

		lappend newData [expr {$pot0 + $pot}]
		incr j
	    }
	}
    }
    
    set grid(data) $newData
    return
}

# Rescale the grid by a constant factor.
proc scalePotential {gridVar factor} {
    upvar $gridVar grid
    
    set newData {}
    foreach pot $grid(data) {
	lappend newData [expr {$factor*$pot}]
    }

    set grid(data) $newData
    return
}
