# potentialGhost.tcl
# Add ghost points.
# Requires gridForce.tcl.
# Author: Jeff Comer <jcomer2@illinois.edu>

if {$argc != 6} {
    puts "This script requires four arguments: "
    puts "  potential map file name (3 columns: x y u)"
    puts "  system size x"
    puts "  system size y"
    puts "  system size z"
    puts "  grid spacing"
    puts "  output file name"
    puts ""
    exit
}

# Parameters:
set inFile [lindex $argv 0]
set lx [lindex $argv 1]
set ly [lindex $argv 2]
set lz [lindex $argv 3]
set dl [lindex $argv 4]
set outFile [lindex $argv end]

source vector.tcl
source gridForce.tcl

# Read an "x y u" lines from a file.
proc readPotential {fileName} {
    set ret {}
    set in [open $fileName r]
    foreach line [split [read $in] \n] {
	if {[string length $line] < 1} {continue}
	if {[string equal [string index $line 0] "\#"]} {continue}

	set tok [concat $line]
	lappend ret [lrange $tok 0 3]
    }
    close $in
    
    return $ret
}

# Return the positions of the 4 closest points.
proc neighboringPoints {rx ry xyuList lx ly} {
    set big 1e40
    set topFour {}
    for {set i 0} {$i < 8} {incr i} {
	lappend topFour [list 0 0 0 $big]
    }

    # Find the nearest four points to rx, ry.
    foreach p $xyuList {
	foreach {px py pu} $p {break}

	# Wrap on the periodic boundaries.
	set dx [expr $rx - $px]
	if {$dx <= -0.5*$lx} {
	    set dx [expr $dx + $lx]
	    set px [expr $px - $lx]
	}
	if {$dx > 0.5*$lx} {
	    set dx [expr $dx - $lx]
	    set px [expr $px + $lx]
	}
	
	set dy [expr $ry - $py]
	if {$dy <= -0.5*$ly} {
	    set dy [expr $dy + $ly]
	    set py [expr $py - $ly]
	}
	if {$dy > 0.5*$ly} {
	    set dy [expr $dy - $ly]
	    set py [expr $py + $ly]
	}
	
	# Compute the distance^2.
	set weight [expr $dx*$dx + $dy*$dy]
	
	# See if this point belongs in the top four.
	if {$weight < [lindex $topFour end end]} {
	    lset topFour end [list $px $py $pu $weight]

	    # Keep the top four sorted.
	    # Still O([llength $xyList])!!!
	    set topFour [lsort -real -decreasing -index end $topFour]
	}
    }

    # Return just the (x,y,u).
    set ret {}
    foreach item $topFour {
	lappend ret [lrange $item 0 2]
    }

    return $ret
}

proc neighboringInterpolate {rx ry xyuList} {
    set xyudList {}
    foreach p $xyuList {
	foreach {px py pu} $p {break}
	set d [expr $px*$px + $py*$py]
	lappend xyudList [list $px $py $pu $d]
    }
    
    set sortList [lsort -real -increasing -index 3 $xyudList]
    # Now we know the origin and the diagonal.
    set p00 [lrange [lindex $sortList 0] 0 2]
    set p11 [lrange [lindex $sortList end] 0 2]

    # Find the other legs.
    set p10 [lrange [lindex $sortList 1] 0 2]
    set p01 [lrange [lindex $sortList 2] 0 2]
    if {[lindex $p10 0] < [lindex $p10 0]} {
	set temp $p10
	set p10 $p01
	set p01 $temp
    }
    
    # Determine the interpolation weights.
    set wx [expr ($rx - [lindex $p00 0])/([lindex $p11 0]-[lindex $p00 0])]
    set wy [expr ($ry - [lindex $p00 1])/([lindex $p11 1]-[lindex $p00 1])]
    
    # Mix along x.
    set g0 [expr [lindex $p10 2]*$wx + (1.0-$wx)*[lindex $p00 2]]
    set g1 [expr [lindex $p11 2]*$wx + (1.0-$wx)*[lindex $p01 2]]
    
    # Mix along lattice vector y.
    return [expr $g1*$wy + (1.0-$wy)*$g0]
}

puts "Reading $inFile"
set data [readPotential $inFile]
puts "Read [llength $data] points."

newGridFit grid $lx $ly $lz $dl

set grid(data) {}
for {set ix 0} {$ix < $grid(nx)} {incr ix} {
    for {set iy 0} {$iy < $grid(ny)} {incr iy} {
	for {set iz 0} {$iz < $grid(nz)} {incr iz} {
	    set r [getPosition grid [list $ix $iy $iz]]
	    set rx [lindex $r 0]
	    set ry [lindex $r 1]
	    set neigh [neighboringPoints $rx $ry $data $lx $ly]
	    
	    lappend grid(data) [neighboringInterpolate $rx $ry $neigh]
	}
    }
}

puts "Generated [llength $grid(data)] potential values."
writeDx grid $outFile


