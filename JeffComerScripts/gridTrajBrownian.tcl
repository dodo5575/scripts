# Author: jcomer2@illinois.edu
source $env(HOME)/scripts/vector.tcl
source $env(HOME)/scripts/gridForce.tcl

proc compute {name moiety structPrefix dcd dcdFreq outDir startFrame} {
    set gridFile phantom_jehanzeb_zero.dx
    set displayPeriod 200
    set stride 1
    set timestep 1.0
    set selText "name POT"
    if {$startFrame < 1} {set startFrame 0}

    # Input:
    set pdb $structPrefix.pdb

    # Load the grid.
    readDx grid $gridFile
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
    mol load pdb $pdb
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
		incr countSum
	    }

	    # Add to velocity map.
	    if {[llength $pos0] != 0} {
		foreach r $pos r0 $pos0 {
		     # Calculate the displacement (angstroms).
		    set d [vecsub $r $r0]
		    
		    # Handle jumps across the periodic boundaries.
		    # Transform to system space.
		    set d [wrapDiff velX $d]
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
    writeDx count $outDir/number_${name}.dx
    writeDx density $outDir/density_${name}.dx
    writeDx velX $outDir/velx_${name}.dx
    writeDx velY $outDir/vely_${name}.dx
    writeDx velZ $outDir/velz_${name}.dx
    writeDx velCount $outDir/velc_${name}.dx
}
