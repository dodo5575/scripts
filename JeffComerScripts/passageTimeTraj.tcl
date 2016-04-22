# jcomer2@uiuc.edu

proc compute {name moiety structPrefix dcd dcdFreq outDir startFrame} {
    set cutZ0 -4
    set cutZ1 4
    set cutDz [expr $cutZ1-$cutZ0]
    set displayPeriod 200
    set stride 1
    set timestep 1.0
    set selText "name POT"
    #set thresh [expr 1.5e-6*$dcdFreq*$timestep]
    if {$startFrame < 1} {set startFrame 1}
    set thresh 0.0

    # Switch the sense of z for negative electric fields.
    if {[string match "*neg*" $name]} {
	set direction -1.0 
    } else {
	set direction 1.0
    }

    # Input:
    set psf $structPrefix.psf
    set pdb $structPrefix.pdb
    set xsc $structPrefix.xsc
    
    # Get the time change between frames in nanoseconds.
    set dt [expr 1.0e-6*$timestep*$dcdFreq*$stride]

    # Read the system size from the xsc file.
    # Note: This only works for lattice vectors along the axes!
    set in [open $xsc r]
    foreach line [split [read $in] "\n"] {
	if {![string match "#*" $line]} {
	    set param [split $line]
	    puts $param
	    set lx [lindex $param 1]
	    set ly [lindex $param 5]
	    set lz [lindex $param 9]
	    break
	}
    }
    puts "NOTE: The system size is $lx $ly $lz.\n"
    close $in

    # Open the output files.
    set out [open "${outDir}/dur_${name}.dat" w]
    #set outAttempt [open "${outDir}/pass_${name}.dat" w]

    # Load the system.
    mol load psf $psf pdb $pdb
    set sel [atomselect top $selText]

    # Loop over the dcd files.
    set nFrames0 0
    foreach dcdFile $dcd {
	# Load the trajectory.
	animate delete all
	mol addfile $dcdFile type dcd step $stride waitfor all
	set nFrames [molinfo top get numframes]
	puts [format "Reading %i frames." $nFrames]

	set indexList [$sel get index]
	set inTime {}
	foreach ind $indexList {
	    lappend inTime -1
	}

	# Move forward, computing
	# current at each step.
	set eventCount 0
	for {set f $startFrame} {$f < $nFrames} {incr f} {
	    # Get the time in nanoseconds for this frame.
	    set t [expr ($nFrames0+$f)*$dt]
	    
	    molinfo top set frame $f
	    set pos [$sel get {x y z}]

	    molinfo top set frame [expr $f-1]
	    set pos0 [$sel get {x y z}]

	    set newInTime {}
	    foreach r $pos r0 $pos0 t0 $inTime ind $indexList {
		set z [expr $direction*[lindex $r 2]]
		set z0 [expr $direction*[lindex $r0 2]]
		
		if {$t0 < 0.0} {
		    # The ion is outside.
		    if {$z0 < $cutZ0 && $z > $cutZ0} {
			# Entering!
			lappend newInTime $t
		    } else {
			# Still outside.
			lappend newInTime -1
		    }
		} else {
		    # The ion is inside.
		    if {$z0 < $cutZ1 && $z > $cutZ1} {
			# Passage event
			set dur [expr $t - $t0]
			# Write it.
			if {$dur > $thresh} { puts $out "$ind $dur" }
			incr eventCount

			lappend newInTime -1
		    } elseif {$z < $cutZ0} {
			# Bounce event
			lappend newInTime -1
		    } else {
			# Still inside
			lappend newInTime $t0
		    }
		}
	    }
	    set inTime $newInTime
	    
	    # Update the display.
	    if {$f % $displayPeriod == 0} {
		puts [format "FRAME %i: %g %i" $f $t $eventCount]
	    }
	}
	set nFrames0 [expr $nFrames+$nFrames0]
    }

    close $out
    mol delete top
}
