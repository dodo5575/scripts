# Author: Jeff Comer <jcomer2@illinois.edu>

proc trajComputeProc {inFileList procName frameTime outFile} {
    set out [open $outFile w]
    set t0 0.0
    set tim 0.0

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
		if {[llength $frame] > 0} {
		    set tim [lindex [concat $line] 1]
		    # Do the analysis.
		    set t [expr {$tim + $t0}]
		    puts $out [$procName $frame $lastFrame $frameTime $t]
		}
		set lastFrame $frame
		set frame {}
		incr frameCount
	    } else {
		# Append the position data
		set tok [concat $line]
		lappend frame $tok
	    }
	}
	close $inp

	# Set the time offset for the start of the next frame.
	set t0 [expr {$t0 + $tim}]
    }
    close $out
}

proc trajComputePassageTime {inFileList frameTime outFile} {
    # Parameters:
    set x0 -475; # Entrance to channel
    set x1 475; # Exit of channel

    set out [open $outFile w]
    set t0 0.0
    set tim 0.0
    set xCut [expr {0.5*($x1+$x0)}]

    set followList {}

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
		if {[llength $frame] == 0} { continue }

		# Get the time.
		set tim [lindex [concat $line] 1]
		set t [expr {$tim + $t0}]
		
		# Remove particles from the list that have exited.
		set followList1 {}
		foreach follow $followList {
		    foreach {followId followTime} $follow { break }

		    set ind [lsearch -index 4 $frame $followId]
		    # Particle got deleted? Shouldn't happen for appropriate x0 and x1.
		    if {$ind < 0 } { 
			puts "Lost particle $followId"
			continue 
		    }
		    
		    foreach {type x y z id} [lindex $frame $ind] { break }
		    if {$x >= $x0 && $x < $x1} {
			# Still inside. Keep it in the follow list.
			lappend followList1 $follow
		    } else {
			if {$x >= $x1} {
			    # The particle has exited the channel.
			    # Get the transit time.
			    set tau [expr {$t - $followTime}]

			    # Write the transit time.
			    puts $out "$id $followTime $t $tau"
			}
		    }
		}
		
		# Set the new followList.
		set followList $followList1
		set nFollow [llength $followList]

		# Add particles to the list to be followed.
		foreach part $frame {
		    foreach {type x y z id} $part { break }
		    # If it's not already in the follow list and has entered the channel,
		    # add it to the follow list.
		    if {[lsearch -index 0 $followList $id] < 0 && $x >= $x0 && $x < $xCut} {
			# Make sure it was outside the channel in the last frame and didn't appear from nowhere.
			# This avoids problems with trajectories that have been running for a while.
			set lastInd [lsearch -index 4 $lastFrame $id]
			if {$lastInd < 0} { continue }
			set lastX [lindex $lastFrame $lastInd 1]

			if {$lastX < $x0} {
			    lappend followList [list $id $t]
			}
		    }
		}

		#puts "FOLLOWING $nFollow -> [llength $followList] particles"

		set lastFrame $frame
		set frame {}
		incr frameCount
	    } else {
		# Append the position data
		set tok [concat $line]
		lappend frame $tok
	    }
	}
	close $inp

	# Set the time offset for the start of the next frame.
	set t0 [expr {$t0 + $tim}]
    }
    close $out
}


# Stores particle binding trajectories with format:
# id chanEnterTime bindTime0 unbindTime0 ... bindTimeN unbindTimeN exitTime passageTime
# Because it is possible to exit without unbinding:
# id chanEnterTime bindTime0 unbindTime0 ... bindTimeN exitTime passageTime
proc trajComputeSurfaceTime {inFileList frameTime outFile} {
    # Parameters:
    set x0 -475; # Entrance to channel
    set x1 475; # Exit of channel

    set y0 -48; # surface location
    set y1 48; # surface location
    set dy 2; # width of surface

    set out [open $outFile w]
    set t0 0.0
    set tim 0.0
    set xCut [expr {0.5*($x1+$x0)}]
    set yCut [expr {0.5*($y1+$y0)}]

    set followList {}

    # Loop through the trajectory files.
    foreach inFile $inFileList {
	set inp [open $inFile r]

	set frameCount 0
	set frame {}
	set lastFrame {}
	# Go through each line of the file.
	while {[gets $inp line] >= 0} {
	    if {[string match "#*" $line]} {continue}
	    if {[string length $line] < 2} {continue}

	    # Find the ends of the frames.
	    if {[string match "END*" $line]} {
		# End of a frame
		if {[llength $frame] == 0} { continue }

		# Get the time.
		set tim [lindex [concat $line] 1]
		set t [expr {$tim + $t0}]
		puts $t
		
		# Remove particles from the list that have exited.
		set followList1 {}
		foreach follow $followList {
		    foreach {followId followTime} $follow { break }

		    set ind [lsearch -index 4 $frame $followId]
		    # Particle got deleted? Shouldn't happen for appropriate x0 and x1.
		    if {$ind < 0 } { 
			puts "Lost particle $followId"
			continue 
		    }
		    
		    foreach {type x y z id} [lindex $frame $ind] { break }
		    if {$x >= $x0 && $x < $x1} {
			# Still inside. Keep it in the follow list.
			# See if it was bound to the surface.
			if {$followBound($id)} {
			    # It was bound.
			    if {($y > $y0+$dy && $y < $yCut) || ($y < $y1-$dy && $y > $yCut)} {
				# The particle has become unbound.
				set followBound($id) 0
				# Append the unbinding time.
				lappend follow $t
			    }
			} else {
			    # It was unbound.
			    if {$y <= $y0 || $y >= $y1} {
				# It has become bound.
				set followBound($id) 1
				# Append the binding time.
				lappend follow $t
			    }
			}

			# Add it back to the follow list.
			lappend followList1 $follow
		    } else {
			if {$x >= $x1} {
			    # The particle has exited the channel.
			    # Get the passage time.
			    set tau [expr {$t - $followTime}]

			    # Append the final time and passage time.
    			    lappend follow $t
			    lappend follow $tau

			    # Write the passage data.
			    puts $out "$follow"
			}
		    }
		}
		
		# Set the new followList.
		set followList $followList1
		set nFollow [llength $followList]

		# Add particles to the list to be followed.
		foreach part $frame {
		    foreach {type x y z id} $part { break }
		    # If it's not already in the follow list and has entered the channel,
		    # add it to the follow list.
		    if {[lsearch -index 0 $followList $id] < 0 && $x >= $x0 && $x < $xCut} {
			# Make sure it was outside the channel in the last frame and didn't appear from nowhere.
			# This avoids problems with trajectories that have been running for a while.
			set lastInd [lsearch -index 4 $lastFrame $id]
			if {$lastInd < 0} { continue }
			set lastX [lindex $lastFrame $lastInd 1]
			if {$lastX >= $x0} { continue }

			if {$y <= $y0 || $y >= $y1} {
			    # The particle is one the surface.
			    set followBound($id) 1
			    # Add both the entrance time and the bound time.
			    lappend followList [list $id $t $t]
			} else {
			    set followBound($id) 0
			    # Add the entrance time.
			    lappend followList [list $id $t]
			}

		    }
		}

		#puts "FOLLOWING $nFollow -> [llength $followList] particles"

		set lastFrame $frame
		set frame {}
		incr frameCount
	    } else {
		# Append the position data
		set tok [concat $line]
		lappend frame $tok
	    }
	}
	close $inp

	# Set the time offset for the start of the next frame.
	set t0 [expr {$t0 + $tim}]
    }
    close $out
}


proc computeNumber {frame lastFrame frameTime t} {
    set x0 -450
    set x1 450

    set n 0
    foreach f $frame {
	foreach {type x y z id} $f { break }
	if {$x >= $x0 && $x < $x1} { incr n }
    }
    
    return [list $t $n]
}


proc computeFlux {frame lastFrame frameTime t} {
    set cutX 450
    set maxJump 100

    set n 0
     foreach f $frame {
	foreach {type x y z id} $f { break }

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

    return [list $t [expr {double($n)/$frameTime}]]
}

proc computeConcCenter {frame lastFrame frameTime t} {
    set x0 -400
    set x1 400
    set y0 -20
    set y1 20
    set lz 100
    
    set vol [expr {double(($x1-$x0)*($y1-$y0)*$lz)}]
    
    set n 0
    foreach f $frame {
	foreach {type x y z id} $f { break }
	if {$x >= $x0 && $x < $x1 && $y >= $y0 && $y < $y1} { incr n }
    }
    
    return [list $t [expr {$n/$vol/6.0221367e-4}]]
}

proc computeConcFront {frame lastFrame frameTime t} {
    set x0 -650
    set x1 -600
    set y0 -20
    set y1 20
    set lz 100
    
    set vol [expr {double(($x1-$x0)*($y1-$y0)*$lz)}]
    
    set n 0
    foreach f $frame {
	foreach {type x y z id} $f { break }
	if {$x >= $x0 && $x < $x1 && $y >= $y0 && $y < $y1} { incr n }
    }
    
    return [list $t [expr {$n/$vol/6.0221367e-4}]]
}

proc computeNumberBound {frame lastFrame frameTime t} {
    set x0 -450
    set x1 450
    set y0 -45
    set y1 45
    
    set n 0
    foreach item $frame {
	foreach {type x y z id} $item { break }
	if {$x >= $x0 && $x < $x1 && ($y < $y0 || $y > $y1)} { incr n }
    }
    
    return [list $t $n]
}


proc computeNumberBound {frame lastFrame frameTime t} {
    set x0 -450
    set x1 450
    set y0 -45
    set y1 45
    
    set n 0
    foreach item $frame {
	foreach {type x y z id} $item { break }
	if {$x >= $x0 && $x < $x1 && ($y < $y0 || $y > $y1)} { incr n }
    }
    
    return [list $t $n]
}


# BROKEN!
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
