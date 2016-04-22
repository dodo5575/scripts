# Author: Jeff Comer <jcomer2@illinois.edu>

proc computeTrajList {inFileList procName frameTime outFile} {
    set out [open $outFile w]
    set t0 0.0

    # Loop through the trajectory files.
    foreach inFile $inFileList {
	set inp [open $inFile r]

	set frameCount 0
	set frame {}
	set lastFrame {}
	while {[gets $inp line] >= 0} {
	    if {[string match "#*" $line]} {continue}
	    if {[string length $line] < 2} {continue}

	    if {[string match "END*" $line]} {
		# End of a frame
		puts $out [$procName $frame $lastFrame $frameTime $t0]
		set lastFrame $frame
		set frame {}
		incr frameCount
	    } else {
		# Get the time.
		set t [lindex $line 1]

		# Append the position data
		set tok [concat $line]
		lappend frame $tok
	    }
	}
	close $inp

	# Set the time offset for the start of the next frame.
	set t0 [expr {$t0 + $t}]
    }
    close $out
}

# t0 is the time offset
proc computeNumber {frame lastFrame frameTime t0} {
    set x0 -450
    set x1 450

    set n 0
    set tim {}
    foreach f $frame {
	foreach {type tim x y z id} $f { break }
	if {$x >= $x0 && $x < $x1} { incr n }
    }
    
    if {[llength $tim] == 0} { return [list $t0 $n] }
    return [list [expr {$tim+$t0}] $n]
}

# t0 is the time offset
proc computeFlux {frame lastFrame frameTime t0} {
    set cutX 450
    set maxJump 100

    set n 0
    set tim {}
    foreach f $frame {
	foreach {type tim x y z id} $f { break }

	# Ignore most particles.
	if {abs($x-$cutX) > $maxJump} { continue }
	
	# Try to find this particle in the other frame.
	set last {}
	foreach part $lastFrame {
	    if {[lindex $part 5] == $id} {
		set last $part
	    }
	}
	if {[llength $last] == 0} { continue }

	set lastX [lindex $last 2]
	if {$lastX <= $cutX && $x > $cutX} {
	    incr n
	} elseif {$lastX >= $cutX && $x < $cutX} {
	    incr n -1
	}
    }

    if {[llength $tim] == 0} { return [list $t0 double($n)/$frameTime] }
    return [list [expr {$tim+$t0}] [expr {double($n)/$frameTime}]]
}


# t0 is the time offset
proc computeConcCenter {frame lastFrame frameTime t0} {
    set x0 -100
    set x1 100
    set y0 -20
    set y1 20
    set lz 100
    
    set vol [expr {($x1-$x0)*($y1-$y0)*$lz}]
    
    set n 0
    foreach f $frame {
	foreach {type tim x y z id} $f { break }
	if {$x >= $x0 && $x < $x1 && $y >= $y0 && $y < $y1} { incr n }
    }
    

    return [list [expr {$tim+$t0}] $n]
}

proc computeConcProfile {inFileList cutTime frameTime outFile} {
    set x0 -625
    set x1 1000
    set y0 -20
    set y1 20
    set lz 100
    set binX 50

    set binN [expr {int(floor(($x1-$x0)/$binX))}]
    set vol [expr {$binX*($y1-$y0)*$lz}] 
    for {set ix 0} {$ix < $binN} {incr ix} {
	set count($ix) 0
    }
    set nFrames 0

    # Loop through the trajectory files.
    foreach inFile $inFileList {
	set inp [open $inFile r]
	
	# Read the file.
	while {[gets $inp line] >= 0} {
	    if {[string match "#*" $line]} {continue}
	    if {[string length $line] < 2} {continue}

	    if {[string match "END*" $line]} {
		# End of a frame
		incr nFrames
	    } else {
		# Get the time.
		set tok [concat $line]
		set t [lindex $tok 1]
		if {$t < $cutTime} { continue }

		set x [lindex $tok 2]
		set y [lindex $tok 3]
		set z [lindex $tok 4]
		
		if {$y >= $y0 && $y < $y1} {
		    set ix [expr {int(floor(($x - $x0)/$binX))}]
		    if {$ix >= 0 && $ix < $binN} {
			incr count($ix)
		    }
		}
	    }
	}
	close $inp
    }

    # Calculate the concentration.
    set out [open $outFile w]
    for {set ix 0} {$ix < $binN} {incr ix} {
	set x [expr {($ix+0.5)*$binX + $x0}]
	set conc [expr {double($count($ix))/$nFrames/$vol/6.0221367e-4}]
	puts $out "$x $conc"
    }

    close $out
}
