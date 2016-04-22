# jcomer2@illinois.edu
if {$argc < 5} {
    puts "$argv0 name regionName structPrefix outputDir dcdFile0 [dcdFile1...]"
    exit
}
set name [lindex $argv 0]
set regionName [lindex $argv 1]
set structPrefix [lindex $argv 2]
set outDir [lindex $argv 3]
set dcdList [lrange $argv 4 end]

source $env(HOME)/scripts/vector.tcl

proc compute {name regionName structPrefix dcdList outDir} {
    set stride 1
    set timestep 1.0
    set dcdFreq 9984

    set z0 75.0
    set rad0 10.0
    set rad1 40.0
    set trajTime 2.0; # time in nanoseconds for each mini-trajectory

    # Define the regions.
    set rad0Sq [expr {$rad0*$rad0}]
    set rad1Sq [expr {$rad1*$rad1}]    
    set region(inTop) "name P and z > $z0 and x^2+y^2<$rad0Sq"
    set region(outTop) "name P and z > $z0 and x^2+y^2>$rad1Sq"
    set region(inBot) "name P and z < $z0 and x^2+y^2<$rad0Sq"
    set region(outBot) "name P and z < $z0 and x^2+y^2>$rad1Sq"
    set region(potLow) "name POT and z < -75"
    set selText $region($regionName)

    # Input:
    set psf $structPrefix.psf
    set pdb $structPrefix.pdb
    set xsc $structPrefix.xsc
    # Output:
    set outFile $outDir/${regionName}_${name}.inst
   
    # Get the time change between frames in nanoseconds.
    set dt [expr 1.0e-6*$timestep*$dcdFreq*$stride]

    # Load the system.
    mol load psf $psf pdb $pdb
    set basisVectors [getCellFromXsc $xsc]
    foreach n {0 1 2} {
	puts "cellBasisVector[expr {$n+1}] [lindex $basisVectors $n]"
    }
    set sys [matTranspose $basisVectors]
    set sysInv [matInvert $sys]
            
    # Loop over the dcd files.
    set out [open $outFile w]
    set nFrames0 0
    set count 0
    foreach dcd $dcdList {
	# Load the trajectory.
	animate delete all
	mol addfile $dcd type dcd step $stride waitfor all
	set nFrames [molinfo top get numframes]
	puts [format "Reading %i frames." $nFrames]

	set cycleFrames [expr {int(floor($trajTime/$dt))}]
	set cycleNum [expr {$nFrames/$cycleFrames}]
	puts "CYCLENUM $cycleNum"
	
	for {set cycle 0} {$cycle < $cycleNum} {incr cycle} {
	    puts "CYCLE $cycle"

	    # First, find particles to follow.
	    set sel [atomselect top $selText]
	    
	    # The first and last frame of the cycle.
	    set f0 [expr {$cycleFrames*$cycle}]
	    set f1 [expr {$f0 + $cycleFrames}]
	    puts "FRAME $f0 TIME [expr {$dt*($f1+$nFrames0)}]"

	    # Trace the trajectories of the selection.
	    set trajList [follow $sel $f0 $f1 $sys $sysInv]
	    set count [expr {$count + [llength $trajList]}]
	    
	    # Write the trajectories.
	    foreach traj $trajList {
		foreach pos $traj {
		    puts $out "$pos"
		}
		puts $out END
	    }

	    $sel delete
	}

	set nFrames0 [expr {$nFrames0 + $nFrames}]
	puts "WROTE $count trajectories"
    }
    
   
    close $out
    mol delete top
}

proc follow {sel frame0 frame1 sys sysInv} {
    set indList [$sel get index]

    # Get the first frame.
    molinfo top set frame $frame0
    set posList [$sel get {x y z}]
    
    # Initialize the trajectory array.
    foreach pos $posList ind $indList {
	set traj($ind) [list $pos]
    }

    # Loop thourgh the frames and grab the trajectories.
    for {set f [expr {$frame0+1}]} {$f < $frame1} {incr f} {
	molinfo top set frame $f

	set posList [$sel get {x y z}]
	foreach pos $posList ind $indList {
	    set pos0 [lindex $traj($ind) end]
	    # Unwrap the trajectory
	    set d [wrapDiff [vecsub $pos $pos0] $sys $sysInv]
	    lappend traj($ind) [vecadd $pos0 $d]
	}
    }

    # Collect the return data.
    set ret {}
    foreach ind $indList {
	lappend ret $traj($ind)
    }
    return $ret
}

proc getCellFromXsc {xsc} {
    # Read the system size from the xsc file.
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
    close $in
    return [list $lx $ly $lz]
}


proc writeTraj {fileName traj} {
    set out [open $fileName w]
    
    foreach item $traj {
	puts $out $item
    }

    close $out

    return
}

proc wrapReal {x l} {
    set l [expr {double($l)}]
    set image [expr {int(floor($x/$l))}]
    set x [expr {$x - $image*$l}]

    return $x
}

proc wrapDiffReal {x l} {
    set l [expr {double($l)}]
    set image [expr {int(floor($x/$l))}]
    set x [expr {$x - $image*$l}]

    if {$x >= 0.5*$l} { set x [expr {$x - $l}] }
    return $x
}

proc wrapDiff {r sys sysInv} {
    set rl [vecTransform $sysInv $r]
    foreach {x y z} $rl { break }
    
    set x [wrapDiffReal $x 1.0]
    set y [wrapDiffReal $y 1.0]
    set z [wrapDiffReal $z 1.0]

    return [vecTransform $sys [list $x $y $z]]
}

compute $name $regionName $structPrefix $dcdList $outDir
exit
