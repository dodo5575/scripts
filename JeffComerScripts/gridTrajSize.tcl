# jcomer2@uiuc.edu
source vector.tcl
source gridForce.tcl

proc compute {name moiety structPrefix dcd dcdFreq outDir startFrame} {
    set len 1.5
    set kind Cl${len}

    set concFactor [expr {29.922*55.523}]
    if {[string match "Cl*" $kind]} {
	set selText "name CLA"
    } else {
	set selText "name POT"
    }

    #set gridFile best_pmf_no.dx
    set displayPeriod 200
    set stride 1
    set timestep 1.0
    if {$startFrame < 1} {set startFrame 0}    

    # Input:
    set psf $structPrefix.psf
    set pdb $structPrefix.pdb
    set xsc $structPrefix.xsc

    # Read the system size from the xsc file.
    set in [open $xsc r]
    foreach line [split [read $in] "\n"] {
	if {![string match "#*" $line]} {
	    set param [split $line]
	    puts $param
	    set ex [lrange $param 1 3]
	    set ey [lrange $param 4 6]
	    set ez [lrange $param 7 9]
	    break
	}
    }
    close $in

    # Load the grid.
    #readDx grid $gridFile
    newGridBox grid $ex $ey $ez [expr double($len)]
    copyGridDim grid count
    copyGridDim grid density
    copyGridDim grid velX
    copyGridDim grid velY
    copyGridDim grid velZ
    copyGridDim grid velCount

    set basis [getSystemCell grid]
    set basisInv [matInvert $basis]
    puts "Grid dimensions: $grid(nx) puts $grid(ny) puts $grid(nz)"
    puts "$grid(delta)"
  
    # Get the time change between frames in nanoseconds.
    set dt [expr 1.0e-6*$timestep*$dcdFreq*$stride]
   
     # Load the system.
    mol load psf $psf pdb $pdb
    set sel [atomselect top $selText]

    # Loop over the dcd files.
    set nFrames0 0
    set countSum 0
    foreach dcdFile $dcd {
	# Load the trajectory.
	animate delete all
	mol addfile $dcdFile type dcd step $stride waitfor all
	set nFrames [molinfo top get numframes]
	puts [format "Reading %i frames." $nFrames]

	# Move forward, computing
	# current at each step.
	set pos0 {}
	for {set f $startFrame} {$f < $nFrames} {incr f} {
	    # Get the time in nanoseconds for this frame.
	    set t [expr ($nFrames0+$f)*$dt]
	    
	    molinfo top set frame $f
	    set pos [$sel get {x y z}]

	    # Add to the density map.
	    foreach r $pos {
		set j [nearestIndex count $r]		
		set c0 [lindex $count(data) $j]
		lset count(data) $j [expr {$c0 + 1.0}]
	    }
	    incr countSum

	    # Add to velocity map.
	    if {[llength $pos0] != 0} {
		foreach r $pos r0 $pos0 {
		     # Calculate the displacement (angstroms).
		    set d [vecsub $r $r0]
		    
		    # Handle jumps across the periodic boundaries.
		    # Transform to system space.
		    set l [vecTransform $basisInv $d]
		    foreach {lx ly lz} $l { break }
		    if {$lx < -0.5} {set lx [expr {$lx + 1.0}]}
		    if {$lx >= 0.5} {set lx [expr {$lx - 1.0}]}
		    if {$ly < -0.5} {set ly [expr {$ly + 1.0}]}
		    if {$ly >= 0.5} {set ly [expr {$ly - 1.0}]}
		    if {$lz < -0.5} {set lz [expr {$lz + 1.0}]}
		    if {$lz >= 0.5} {set lz [expr {$lz - 1.0}]}
		    set l [list $lx $ly $lz]
		    # Transform back to world space.
		    set d [vecTransform $basis $l]
		    foreach {vx vy vz} $d { break }

		    # Find the midpoint respecting the periodic boundaries.
		    set m [vecadd $r0 [vecscale 0.5 $d]]
		    # Get the data for the corresponding node.
		    set j [nearestIndex velX $m]
		    set vx0 [lindex $velX(data) $j]
		    set vy0 [lindex $velY(data) $j]
		    set vz0 [lindex $velZ(data) $j]
		    set c0 [lindex $velCount(data) $j]
		   
		    # Add this contribution to the node.
		    lset velX(data) $j [expr {$vx0 + $vx}]
		    lset velY(data) $j [expr {$vy0 + $vy}]
		    lset velZ(data) $j [expr {$vz0 + $vz}]
		    lset velCount(data) $j [expr {$c0 + 1.0}]
		}
	    }

	    # Update the display.
	    if {$f % $displayPeriod == 0} {
		puts [format "FRAME %i" $f ]
	    }
	    set pos0 $pos
	}
	set nFrames0 [expr $nFrames+$nFrames0]
    }
    mol delete top
    $sel delete

    # Compute the mean number density.
    set volume [getVolume count]
    set meanDensity {}
    set meanCount {}
    foreach p $count(data) {
	lappend meanCount [expr {$p/$countSum}]
	lappend meanDensity [expr {$p/$countSum/$volume}]
    }
    set count(data) $meanCount
    set density(data) $meanDensity

    # Compute the mean velocities.
    set meanVx {}
    set meanVy {}
    set meanVz {}
    foreach vx $velX(data) vy $velY(data) vz $velZ(data) c $velCount(data) {
	if {$c > 0} {
	    lappend meanVx [expr {$vx/$c/$dt}]
	    lappend meanVy [expr {$vy/$c/$dt}]
	    lappend meanVz [expr {$vz/$c/$dt}]
	} else {
	    lappend meanVx 0.0
	    lappend meanVy 0.0
	    lappend meanVz 0.0
	}
    }
    set velX(data) $meanVx
    set velY(data) $meanVy
    set velZ(data) $meanVz

    # Write the results.
    writeDx count $outDir/number${kind}_${name}.dx
    writeDx density $outDir/density${kind}_${name}.dx
    writeDx velX $outDir/velx${kind}_${name}.dx
    writeDx velY $outDir/vely${kind}_${name}.dx
    writeDx velZ $outDir/velz${kind}_${name}.dx
    writeDx velCount $outDir/velc${kind}_${name}.dx
}
