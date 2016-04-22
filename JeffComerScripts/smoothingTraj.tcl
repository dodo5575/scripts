# jcomer2@uiuc.edu

source vector.tcl

proc makeCell {cellVar a b c} {
    upvar $cellVar cell
    set cell(basis) [matTranspose [list $a $b $c]]
    set cell(basisInv) [matInvert $cell(basis)]
    set cell(origin) [vecscale -0.5 [vecAdd $a [vecAdd $b $c]]]
    return
}

# Cell has members cell(basis), cell(basisInv), cell(origin)
# cell(basis) is a 3*3 matrix of the cell basis vectors as columns.
#set cell(basisInv) [matInvert $cell(basis)]
proc wrap {cellVar r} {
    upvar $cellVar cell

    set l [vecTransform $cell(basisInv) [vecsub $r $cell(origin)]]
    foreach {x y z} $l { break }
    
    set lwx [expr {$x - int(floor($x))}]
    set lwy [expr {$y - int(floor($y))}]
    set lwz [expr {$z - int(floor($z))}]
    set rw [vecAdd [vecTransform $cell(basis) [list $lwx $lwy $lwz]]  $cell(origin)]
    return $rw 
}

proc wrapDisplace {cellVar d} {
    upvar $cellVar cell

    set l [vecTransform $cell(basisInv) $d]
    foreach {x y z} $l { break }
    
    set lwx [expr {$x - floor($x)}]
    if {$lwx >= 0.5} { set lwx [expr {$lwx-1.0}] }

    set lwy [expr {$y - floor($y)}]
    if {$lwy >= 0.5} { set lwy [expr {$lwy-1.0}] }

    set lwz [expr {$z - floor($z)}]
    if {$lwz >= 0.5} { set lwz [expr {$lwz-1.0}] }

    set rw [vecTransform $cell(basis) [list $lwx $lwy $lwz]]
    return $rw 
}

proc wrapList {cellVar posList} {
    upvar $cellVar cell
    
    set retList {}
    foreach pos $posList {
	lappend retList [wrap cell $pos]
    }
    return $retList
}

proc unwrapTrajectory {cellVar posListTrajVar} {
    upvar $cellVar cell
    upvar $posListTrajVar posListTraj

    set newListTraj {}
    set posList0 [lindex $posListTraj 0]
    foreach posList $posListTraj {
	set newList {}
	foreach pos $posList pos0 $posList0 {
	    set d [wrapDisplace cell [vecSub $pos $pos0]]
	    lappend newList [vecAdd $pos0 $d]
	}
	lappend newListTraj $newList
	set posList0 $newList
    }
    return $newListTraj
}

# Assumes that the position list is unwrapped!
proc meanPosition {posListList} {
    set sumList {}
    foreach pos [lindex $posListList 0] {
	lappend sumList [vecZero]
    }

    # Add together all the time steps.
    foreach posList $posListList {
	set newSumList {}
	foreach pos $posList sum $sumList {
	    lappend newSumList [vecAdd $pos $sum]
	}
	set sumList $newSumList
    }

    # Divide by the number of time steps.
    # Wrap the resulting positions.
    set n [llength $posListList]
    set avgList {}
    foreach sum $sumList {
	lappend avgList [vecscale [expr {1.0/$n}] $sum]
    }
    return $avgList
}

proc compute {name moiety structPrefix dcdList dcdFreq outDir stride windowRad} {
    set timestep 1.0
    set selText all
    set outFile smooth${windowRad}_${name}.dcd

    # Input:
    set psf $structPrefix.psf
    set pdb $structPrefix.pdb
    set xsc $structPrefix.xsc

    # Read the system size from the xsc file.
    # Note: This only works for lattice vectors along the axes!
    set in [open $xsc r]
    foreach line [split [read $in] "\n"] {
	if {![string match "#*" $line]} {
	    set param [split $line]
	    puts $param
	    set lx [lrange $param 1 3]
	    set ly [lrange $param 4 6]
	    set lz [lrange $param 7 9]
	    break
	}
    }
    puts "lx: $lx"
    puts "ly: $ly"
    puts "lz: $lz"
    close $in
    
    # Make the unit cell.
    makeCell cell $lx $ly $lz

    # Load the system.
    mol load psf $psf pdb $pdb
    set sel [atomselect top $selText]
    set all [atomselect top all]
    puts "Loaded the structure `$psf' `$pdb'."
    animate delete all
 
    # Loop over the dcd files.
    set nFrames0 0
    foreach dcdFile $dcdList {
	# Load the trajectory.
	mol addfile $dcdFile type dcd step $stride waitfor all
	set nFrames [molinfo top get numframes]
	puts [format "Reading %i frames." [expr {$nFrames-$nFrames0}]]
	incr nFrames0 $nFrames
    }
    set totalFrames [molinfo top get numframes]
    puts "Loaded [llength $dcdList] dcd files."
    
    # Make a list of all atom positions for all frames.
    set posListTraj {}
    for {set f 0} {$f < $totalFrames} {incr f} {
	molinfo top set frame $f
	lappend posListTraj [$sel get {x y z}]
    }
    puts "Extracted the trajectory."

    # Unwrap the trajectory.
    set unwrapListTraj [unwrapTrajectory cell posListTraj]
    puts "Unwrapped the trajectory."
    
    # Regenerate each frame.
    for {set f 0} {$f < $totalFrames} {incr f} {
	molinfo top set frame $f
	
	# Find the indices over which to average.
	set f0 [expr {$f - $windowRad}]
	if {$f0 < 0} { set f0 0 }
	set f1 [expr {$f + $windowRad}]
	if {$f1 >= $totalFrames} { set f1 [expr {$totalFrames-1}] }
	
	# Get the mean positions,
	# wrapping to the point nearest the original point.
	set buffer [lrange $unwrapListTraj $f0 $f1]
	set now [lindex $posListTraj $f]
	set meanPosList [meanPosition $buffer]
	# Wrap it back to the original cell.
	set placementTraj [list $now $meanPosList]
	set finalTraj [unwrapTrajectory cell placementTraj]
	$sel set {x y z} [lindex $finalTraj 1]
	#$sel set {x y z} [lindex $unwrapListTraj $f]
	
	puts "Smoothed frame $f."
    }

    animate write dcd $outDir/$outFile beg 0 end -1 sel $all top
    puts "Wrote the output dcd."

    $sel delete
    $all delete
    mol delete top
}
