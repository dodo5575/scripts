# jcomer2@uiuc.edu
proc writeTraj {fileName traj} {
    set out [open $fileName w]
    
    foreach item $traj {
	puts $out $item
    }

    close $out

    return
}

proc wrap {d l} {
    foreach {dx dy dz} $d { break }
    foreach {lx ly lz} $l { break }
    if {$dx < -0.5*$lx} { set dx [expr $dx + $lx] }
    if {$dx > 0.5*$lx} { set dx [expr $dx - $lx] }
    if {$dy < -0.5*$ly} { set dy [expr $dy + $ly] }
    if {$dy > 0.5*$ly} { set dy [expr $dy - $ly] }
    if {$dz < -0.5*$lz} { set dz [expr $dz + $lz] }
    if {$dz > 0.5*$lz} { set dz [expr $dz - $lz] }
    
    return [list $dx $dy $dz]
}

proc compute {name moiety structPrefix dcdList dcdFreq outDir startFrame} {
    set displayPeriod 200
    set stride 1
    set timestep 1.0
    set lx 40
    set ly 40
    set lz 72
    set limitZ 30
    set selText "name POT and not resid 1"
    # Input:
    set psf $structPrefix.psf
    set pdb $structPrefix.pdb
    # Output:
    set outPrefix $outDir/follow_${name}

    # Get the time change between frames in nanoseconds.
    set dt [expr 1.0e-6*$timestep*$dcdFreq*$stride]

    # Load the system.
    mol load psf $psf pdb $pdb
    set sel [atomselect top $selText]
    set sysSize [list $lx $ly $lz]

    # Loop over the dcd files.
    set nFrames0 0
    set count 0
    foreach dcd $dcdList {
	# Load the trajectory.
	animate delete all
	mol addfile $dcd type dcd step $stride waitfor all
	set nFrames [molinfo top get numframes]
	puts [format "Reading %i frames." $nFrames]

	set followIon -1
	set pos0 [veczero]
	set followTraj {}

	# Move forward computing at every step.
	for {set f 0} {$f < $nFrames} {incr f} {
	    molinfo top set frame $f
	    set t [expr ($nFrames0+$f)*$dt]
	    set posList [$sel get {x y z}]

	    if {$followIon > 0} {
		set r [lindex $posList $followIon]
		if {abs([lindex $r 2]) < $limitZ} {
		    # Write the pending trajectory.
		    writeTraj ${outPrefix}${count}.dat $followTraj
		    incr count
		    
		    # Get ready for a new trajectory.
		    set followIon -1
		    set followTraj {}
		} else {
		    # Store this position in the trajectory.
		    set r0 [lrange [lindex $followTraj end] 1 3]
		    # Unwrap the trajectory.
		    set d [wrap [vecsub $r $r0] $sysSize]
		    set r1 [vecadd $r0 $d]
		    lappend followTraj [concat $t $r1]
		}
	    } else {
		# Find the best choice for a ion to follow.
		set i 0
		set maxZ [expr abs([lindex $posList 0 2])]
		set maxIndex 0
		foreach pos $posList {
		    set z [expr abs([lindex $pos 2])]
		    if {$z > $maxZ} {
			set maxZ $z 
			set maxIndex $i
		    }
		    incr i
		}
		
		# Should we follow this ion?
		if {$maxZ > $limitZ} {
		    # Follow it.
		    set followIon $maxIndex
		    set followTraj [list [concat $t [lindex $posList $maxIndex]]]
		}
	    }

	    # Write the time and distance.
	    if {$f % $displayPeriod == 0} {
		puts [format "FRAME %i: %i events" $f $count]
	    }
	}
	set nFrames0 [expr $nFrames+$nFrames0]
    }

    $sel delete
    mol delete top
}
