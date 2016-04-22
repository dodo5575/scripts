# jcomer2@uiuc.edu
 
# Return the current and number of charge carriers in the selection for the frame.
# Does not compensate for jumps over the periodic boundary, although
# they should be ignored for an appropriate choice of cutZ0 and cutZ1.
proc computeCurrent {frameCurr sel charge ly dt cutZ0 cutZ1} {
    set frameLast [expr {$frameCurr-1}]
    
    # Get the current position of the selection.
    molinfo top set frame $frameCurr
    set z1List [$sel get z]
    set y1List [$sel get y]
    
    # Get the last position of the selection.
    molinfo top set frame $frameLast
    set z0List [$sel get z]
    set y0List [$sel get y]
    
    # Get the number of charge carriers.
    set num [$sel num]
    if {$num == 0} { return [list 0 0.0]}
    
    # Find the displacements in the z-direction and compute the current.
    set currentY 0.0
    foreach z0 $z0List z1 $z1List y0 $y0List y1 $y1List {
	# Ignore carriers outside the range for both times.
	if {($z0 < $cutZ0 || $z0 > $cutZ1) && ($z1 < $cutZ0 || $z1 > $cutZ1)} {
	    continue
	}

	set az0 $z0
	set az1 $z1
	if {$z0 < $cutZ0} { set az0 $cutZ0 }
	if {$z1 < $cutZ0} { set az1 $cutZ0 }
	if {$z0 > $cutZ1} { set az0 $cutZ1 }
	if {$z1 > $cutZ1} { set az1 $cutZ1 }
	
	if {$z1-$z0 == 0.0} { continue }
	
	set m [expr {($y1-$y0)/($z1-$z0)}]
	set ay0 [expr {$m*($az0-$z0) + $y0}]
	set ay1 [expr {$m*($az1-$z0) + $y0}]
	set dy [expr {$ay1 - $ay0}]

	# Compute the current in nanoamperes.
	set currentY [expr {$currentY + $charge*$dy/($ly*$dt)*1.60217733e-1}]
    }
    
    return [list $num $currentY]
}

proc compute {name moiety structPrefix dcd dcdFreq outDir stride} {
    set lenY [expr {110.0*sqrt(3.0)/2.0}]
    set cutZ0 48.0
    set cutZ1 78.0
    set displayPeriod 200
    set startFrame 1
    set timestep 1.0
    set selText [list "name POT" "name CLA"]
    set charge [list 1.0 -1.0]
    set nameList [list "K" "Cl"]
    #set startFrame 600
    if {$stride < 1} { set stride 1 }

    # Input:
    set psf $structPrefix.psf
    set pdb $structPrefix.pdb
    #set xsc $structPrefix.xsc

    
    # Get the time change between frames in nanoseconds.
    set dt [expr {1.0e-6*$timestep*$dcdFreq*$stride}]

    # Open the output files.
    set out {}
    foreach n $nameList {
	lappend out [open "${outDir}/curr${n}_${name}.dat" w]
    }
    set outTotal [open "${outDir}/curr_${name}.dat" w]

    # Load the system.
    mol load psf $psf pdb $pdb
    set sel {}
    foreach st $selText {
	lappend sel [atomselect top $st]
    }

    # Loop over the dcd files.
    set nFrames0 0
    foreach dcdFile $dcd {
	# Load the trajectory.
	animate delete all
	mol addfile $dcdFile type dcd step $stride waitfor all
	set nFrames [molinfo top get numframes]
	puts [format "Reading %i frames." $nFrames]

	# Move forward, computing
	# current at each step.
	for {set f $startFrame} {$f < $nFrames} {incr f} {
	    # Get the time in nanoseconds for this frame.
	    set t [expr {($nFrames0+$f+0.5)*$dt}]

	    # Compute the current for each selection.		
	    set currentZTotal 0.0
	    foreach s $sel q $charge o $out {
		#Write the number of carriers and the current for this selection.
		set data [computeCurrent $f $s $q $lenY $dt $cutZ0 $cutZ1]
		set currentZ [lindex $data 1]
		puts $o "$t $currentZ"
		
		set currentZTotal [expr {$currentZTotal + $currentZ}]
	    }
	    
	    # Write the total current.
	    puts $outTotal "$t $currentZTotal"
	    
	    # Update the display.
	    if {$f % $displayPeriod == 0} {
		puts -nonewline [format "FRAME %i: " $f]
		puts "$t $currentZTotal"
	    }
	}
	set nFrames0 [expr {$nFrames+$nFrames0}]
    }

    foreach s $sel o $out {
	$s delete
	close $o
    }
    close $outTotal
    mol delete top
}
