# Author: Jeff Comer <jcomer2@illinois.edu>

proc countNeighbors {sysVar ix iy nx ny} {
    upvar $sysVar sys

    # Wrap the neighbors.
    set ix0 [expr $ix - 1]
    if {$ix0 < 0} { set ix0 [expr $nx-1] }
    set ix1 [expr $ix + 1]
    if {$ix1 >= $nx} { set ix1 0 }
    set iy0 [expr $iy - 1]
    if {$iy0 < 0} { set iy0 [expr $ny-1] }
    set iy1 [expr $iy + 1]
    if {$iy1 >= $ny} { set iy1 0 }

    # Count the eight neighbors.
    set count 0
    if {$sys($ix0,$iy0)} { incr count }
    if {$sys($ix,$iy0)} { incr count }
    if {$sys($ix1,$iy0)} { incr count }
    if {$sys($ix0,$iy1)} { incr count }
    if {$sys($ix,$iy1)} { incr count }
    if {$sys($ix1,$iy1)} { incr count }
    if {$sys($ix0,$iy)} { incr count }
    if {$sys($ix1,$iy)} { incr count }

    return $count
}

proc wrap {ix iy nx ny} {
    if {$ix < 0} { set ix [expr $ix + $nx] }
    if {$ix >= $nx} { set ix [expr $ix - $nx] }
    if {$iy < 0} { set iy [expr $iy + $ny] }
    if {$iy >= $ny} { set iy [expr $iy - $ny] }
}

proc loadPipes {sysVar fileName} {
    upvar $sysVar sys 
    
    set in [open $fileName r]
    set iy 0
    while {[gets $in line] >= 0} {
	if {[string length $line] <= 1} { continue }
        if {[string match "#*" $line]} { continue }

	set tok [concat $line]
	set nx [llength $tok]

	for {set ix 0} {$ix < $nx} {incr ix} {
	    set sys($ix,$iy) [lindex $tok $ix]
	}

	incr iy
    }
    set ny $iy
    close $in

    return [list $nx $ny]
}

# Nothing 0
# Other 1
# End 2
# Horizontal 3
# Vertical 4
proc preparePipes {sys nx ny} {
    upvar $sysVar sys

    for {set ix 0} {$ix < $nx} {incr ix} {
	for {set iy 0} {$iy < $ny} {incr iy} {
	    set c [countNeighbors sys $ix $iy $nx $ny]
	    if {$c == 1} {
		# End.
		set set sys($ix,$iy) 2
	    } else {
		set x0 [wrap [expr $ix-1] $nx $ny]
		set x1 [wrap [expr $ix+1] $nx $ny]
		set y0 [wrap [expr $iy-1] $nx $ny]
		set y1 [wrap [expr $iy+1] $nx $ny]

		if {$sys($x0,$iy)>0 && $sys($x1,$iy)>0} {
		    # Horizontal.
		    set sys($ix,$iy) 3
		} elseif {$sys($ix,$y0)>0 && $sys($ix,$y1)>0} {
		    # Vertical.
		    set sys($ix,$iy) 4
		}
	    }
	    
	}
    }
    return
}

proc drawPipes {sys nx ny origin size} {
    graphics top delete all
    set sphereRad [expr 0.6*$size]
    set cylRad [expr $0.2*$size]

    for {set ix 0} {$ix < $nx} {incr ix} {
	for {set iy 0} {$iy < $ny} {incr iy} {
	    set x [expr $ix*$size]
	    set y [expr $iy*$size]
	    set pos [vecadd [list $x $y 0] $origin]

	    set state $sys($ix,$iy)
	    
	    if {$state == 0} { continue }
	    if {$state == 1} {
		# Other
		graphics top sphere $pos radius $sphereRad resolution 12
	    } elseif {$state == 2} {
		# End
	    } elseif {$state == 3} {
		# Horizontal
		set pos0 [vecsub $pos [list [expr -$size] 0 0]]
		set pos1 [vecsub $pos [list $size 0 0]]
		cylinder $pos0 $pos1 radius $cylRad resolution 12
	    } elseif {$state == 3} {
		# Vertical
		set pos0 [vecsub $pos [list 0 [expr -$size] 0]]
		set pos1 [vecsub $pos [list 0 $size 0]]
		cylinder $pos0 $pos1 radius $cylRad resolution 12
	    }
	}
    }
}
