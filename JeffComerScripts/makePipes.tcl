# Author: Jeff Comer <jcomer2@illinois.edu>
set nx 100
set ny 50
set n [expr $nx*$ny]
set steps 100
set start [expr 2*$nx/3]
set turnProb 0.1
set outPrefix big_maze

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

# Make a cellular automaton step.
proc stepMaze {sysVar nx ny} {
    upvar $sysVar sys

    # Count the neighbors.
    for {set ix 0} {$ix < $nx} {incr ix} {
	for {set iy 0} {$iy < $ny} {incr iy} {
	    set count($ix,$iy) [countNeighbors sys $ix $iy $nx $ny]
	}
    }

    for {set ix 0} {$ix < $nx} {incr ix} {
	for {set iy 0} {$iy < $ny} {incr iy} {
	    set c $count($ix,$iy)
	    
	    if {$sys($ix,$iy) == 0} {
		# Check for birth.
	        if {$c == 3} { set sys($ix,$iy) 1 }
	    } else {
		# Check for death.
		if {$c < 1 || $c > 5} { set sys($ix,$iy) 0 }
	    }
	}
    }

    return
}

proc writeSystem {fileName sysVar nx ny} {
    upvar $sysVar sys

    set out [open $fileName w]
    for {set iy 0} {$iy < $ny} {incr iy} { 
	for {set ix 0} {$ix < $nx} {incr ix} {
	    puts -nonewline $out $sys($ix,$iy)
	}
	puts $out ""
    }
    close $out
}

proc clear {sysVar nx ny} {
    upvar $sysVar sys
    
    for {set ix 0} {$ix < $nx} {incr ix} {
	for {set iy 0} {$iy < $ny} {incr iy} {	
	    set sys($ix,$iy) 0
	}
    }
    return
}

proc initial {sysVar nx ny start turnProb} {
    upvar $sysVar sys
    set x0 [expr $nx/2 - $start/2]
    set x1 [expr $x0 + $start]
    set iy [expr $ny/2]
    for {set ix $x0} {$ix < $x1} {incr ix} {
	set sys($ix,$iy) 1

	# Decide whether to turn or not.
	set r [expr rand()]
	if {$r < 0.5*$turnProb} {
	    incr iy -1
	    if {$iy < 0} { set iy [expr $ny-1] }
	    set sys($ix,$iy) 1
	} elseif {$r < $turnProb} {
	    incr iy
	    if {$iy >= $ny} { set iy 0 }
	    set sys($ix,$iy) 1
	}
    }
    return
}


# Set the initial condition.
clear pipe $nx $ny
set cx [expr $nx/2]
set cy [expr $ny/2]
set pipe($cx,$cy) 1
set pipe([expr $cx-1],$cy) 1
set pipe([expr $cx+1],$cy) 1
initial pipe $nx $ny $start $turnProb

# Run Maze.
writeSystem ${outPrefix}0.cell pipe $nx $ny
for {set s 0} {$s < $steps} {incr s} {
    stepMaze pipe $nx $ny
    if {$s % 20 == 0} { puts "step $s" }
}
writeSystem ${outPrefix}.cell pipe $nx $ny
