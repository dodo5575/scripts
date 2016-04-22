# computeGridForce.tcl
# Grid force procedures
# Author: Jeff Comer <jcomer2@illinois.edu>

# Returns the potential value at each point in gridVar(data),
# the number of grid points along each
# cell axis in gridVar(nx,ny,nz), the lattice vectors in
# gridVar(delta) (3x3 matrix), and the cell origin in
# gridVar(origin) (3-vector).
# Data order is z fast, y medium, and x slow.
proc readDx {gridVar fileName} {
    upvar grid $gridVar
    set in [open $fileName r]

    set items 0
    set n 0
    set grid(data) {}
    set grid(delta) {}

    foreach line [split [read $in] "\n"] {
	set tok [concat $line]
	
	# Ignore comments and blank lines.
	if {[string match "\#*" $tok]} {continue}
	if {[llength $line] < 1} {continue}

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
}

# Writes a dx file given a grid structure as returned by readDx.
proc writeDx {gridVar fileName} {
    upvar grid $gridVar
    
    foreach {ox oy oz} $grid(origin) {break}
    set delta [join $grid(delta)]
    foreach {exx exy exz eyx eyy eyz ezx ezy ezz} $delta {break}

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




