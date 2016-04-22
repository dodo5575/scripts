# jcomer2@uiuc.edu

# Return the current and number of charge carriers in the selection for the frame.
# Does not compensate for jumps over the periodic boundary, although
# they should be ignored for an appropriate choice of cutZ0 and cutZ1.
# cyl has cyl(s0), cyl(ds), cyl(ns). Also, cyl(0), cyl(1)...corresponding to radial bins.
# Make sure you zero cyl(0..n) ahead of time.
proc computeCurrent {cylGridVar frameCurr sel charge dt cutZ0 cutZ1} {
    upvar $cylGridVar cyl
    set frameLast [expr {$frameCurr-1}]
    set lz [expr {$cutZ1-$cutZ0}]

    # Get the current position of the selection.
    molinfo top set frame $frameCurr
    set pos1 [$sel get {x y z}]

    # Get the last position of the selection.
    molinfo top set frame $frameLast
    set pos0 [$sel get {x y z}]   

    # Find the displacements in the z-direction and compute the current (pA).
    foreach r0 $pos0 r1 $pos1 {
	# Get the z coordinates.
	set z0 [lindex $r0 2]
	set z1 [lindex $r1 2]

	# Ignore carriers outside the range for both times.
	if {($z0 < $cutZ0 || $z0 > $cutZ1) && ($z1 < $cutZ0 || $z1 > $cutZ1)} {
	    continue
	}

	# Get the radial coordinate.
	set mid [vecscale 0.5 [vecadd $r0 $r1]]
	foreach {x y z} $mid { break }
	set s [expr {sqrt($x*$x +  $y*$y)}]
	set x [expr {($s-$cyl(s0))/$cyl(ds)}]
	set is [expr {int(floor($x))}]
	
	# Cut pieces outside of the range.
	if {$z0 < $cutZ0} { set z0 $cutZ0 }
	if {$z1 < $cutZ0} { set z1 $cutZ0 }
	if {$z0 > $cutZ1} { set z0 $cutZ1 }
	if {$z1 > $cutZ1} { set z1 $cutZ1 }
	set dz [expr {$z1-$z0}]


	# Compute the current in picoamperes.
	set curr [expr {$charge*$dz/($lz*$dt)*1.60217733e2}]
	set cyl($is) [expr {$cyl($is) + $curr}]
	set cyl(s,$is) [expr {$cyl(s,$is) + $curr*$curr}]
    }
}

proc compute {name moiety structPrefix dcdList dcdFreq outDir startFrame} {
    set displayPeriod 200
    set timestep 1.0
    set stride 1

    set posName ionK
    set posText "name POT"
    set posCharge 1.0
    
    set negName ionCl
    set negText "name CLA"
    set negCharge -1.0

    if {$stride <= 0} {set stride 1}
    if {$startFrame <= 0} {set startFrame 1}
    set avogadro 6.0221367e23

    # The geometry:
    set z0 -12
    set z1 12
    set s0 0.0
    set s1 30.0
    set ns 60

    # Load the geometry variables.
    set pos(s0) $s0
    set pos(ds) [expr {($s1-$s0)/$ns}]
    set pos(ns) $ns
    set neg(s0) $s0
    set neg(ds) [expr {($s1-$s0)/$ns}]
    set neg(ns) $ns

    # Zero the sums.
    for {set is 0} {$is < $ns} {incr is} {
	set pos($is) 0.0
	set neg($is) 0.0
	set pos(s,$is) 0.0
	set neg(s,$is) 0.0
    }

    # Input:
    set pdb $structPrefix.pdb
  
    # Get the time change between frames in nanoseconds.
    set dt [expr {1.0e-6*$timestep*$dcdFreq*$stride}]

    # Load the system.
    mol load pdb $pdb
    set posSel [atomselect top $posText]
    set negSel [atomselect top $negText]

    # Loop over the dcd files.
    set nFrames0 0
    set totalCount 0
    foreach dcd $dcdList {
	# Load the trajectory.
	animate delete all
	mol addfile $dcd type dcd step $stride waitfor all
	set nFrames [molinfo top get numframes]
	puts [format "Reading %i frames." $nFrames]

	# Move forward computing at every step.
	for {set f $startFrame} {$f < $nFrames} {incr f} {
	    incr totalCount
	    
	    computeCurrent pos $f $posSel $posCharge $dt $z0 $z1
	    computeCurrent neg $f $negSel $negCharge $dt $z0 $z1

	    # Write the time and distance.
	    if {$f % $displayPeriod == 0} {
		puts -nonewline [format "FRAME %i: " $f]
		puts "$totalCount"
	    }
	}
	set nFrames0 [expr $nFrames+$nFrames0]
    }

    # Write the results.
    set posOut [open $outDir/${posName}_${name}.dat w]
    set negOut [open $outDir/${negName}_${name}.dat w]
    for {set is 0} {$is < $ns} {incr is} {
	
	# Get the count.
	set n [expr {double($totalCount)}]
	
	if {$n <= 0} {
	    set posCurr 0.0
	    set posErr 0.0
	    set negCurr 0.0
	    set negErr 0.0
	} else {
	    set posCurr [expr {$pos($is)/$n}]
	    set posErr [expr {sqrt(($pos(s,$is) - $pos($is)*$pos($is)/$n)/($n-1)/$n)}]
	    set negCurr [expr {$neg($is)/$n}]
	    set negErr [expr {sqrt(($neg(s,$is) - $neg($is)*$neg($is)/$n)/($n-1)/$n)}]
	}

	set sm [expr {$s0 + $pos(ds)*($is+0.5)}]
	puts $posOut "$sm $posCurr $posErr"
	puts $negOut "$sm $negCurr $negErr"
    }

    close $posOut
    close $negOut
    $posSel delete
    $negSel delete
    mol delete top
}
