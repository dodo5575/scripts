# jcomer2@illinois.edu
if {$argc < 8} {
    puts "$argv0 name regionName trajTime structPrefix outputDir dcdFreq stride dcdFile0 [dcdFile1...]"
    exit
}
set name [lindex $argv 0]
set regionName [lindex $argv 1]
set trajTime [lindex $argv 2]; # time in ns for each mini-trajectory
set structPrefix [lindex $argv 3]
set outDir [lindex $argv 4]
set dcdFreq [lindex $argv 5]
set stride [lindex $argv 6]
set dcdList [lrange $argv 7 end]

source $env(HOME)/scripts/vector.tcl

proc compute {name regionName trajTime structPrefix dcdFreq dcdList outDir stride} {
    if {$stride < 1} { set stride 1 }
    set timestep 1.0
    set xsc ../extra_layer_dopc_neg.restart.xsc

    set zBelow -65
    # Define the regions.
    set region(potFree) "name POT and z < $zBelow"
    set region(chlFree) "name CLA and z < $zBelow"
    set selText $region($regionName)

      # Input:
    set psf $structPrefix.psf
    set pdb $structPrefix.pdb
    #set xsc $structPrefix.xsc
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
    animate delete all
    foreach dcd $dcdList {
	# Load the trajectory.
	mol addfile $dcd type dcd step $stride waitfor all
    }
    set nFrames [molinfo top get numframes]
    puts [format "Reading %i frames." $nFrames]

    set cycleFrames [expr {int(floor($trajTime/$dt))}]
    set cycleNum [expr {$nFrames/$cycleFrames}]
    puts "CYCLENUM $cycleNum"

    for {set cycle 0} {$cycle < $cycleNum} {incr cycle} {
	puts "CYCLE $cycle"
	# The first and last frame of the cycle.
	set f0 [expr {$cycleFrames*$cycle}]
	set f1 [expr {$f0 + $cycleFrames}]

	# Find particles to follow.
	molinfo top set frame $f0
	set sel [atomselect top $selText]
	set resList [lsort -unique [$sel get {segname resid}]]
	$sel delete

	set tim [expr {$dt*($f1+$nFrames0)}]
	puts "FRAME $f0 TIME $tim NUM [llength $resList]"

	# Trace the trajectories of each residue.
	set trajList {}
	foreach res $resList {
	    foreach {segName resId} $res { break }
	    set s [atomselect top "resid $resId and segname $segName"]
	    lappend trajList [followCenter $s $f0 $f1 $sys $sysInv]
	    $s delete
	}
	set count [expr {$count + [llength $trajList]}]
	
	# Write the trajectories.
	foreach traj $trajList {
	    foreach pos $traj {
		puts $out "$pos"
	    }
	    puts $out "END $tim"
	}

	puts "WROTE $count trajectories"
    }

    close $out
    mol delete top
}

proc followCenter {sel frame0 frame1 sys sysInv} {
    # Get the first position.
    molinfo top set frame $frame0
    set pos0 [measure center $sel weight mass]
    set traj [list $pos0]

    # Loop thourgh the frames and grab the trajectories.
    for {set f [expr {$frame0+1}]} {$f < $frame1} {incr f} {
	# Get the center of mass at this frame.
	molinfo top set frame $f
	set r [measure center $sel weight mass]

	# Unwrap the trajectory
	set d [wrapDiff [vecsub $r $pos0] $sys $sysInv]
	set pos [vecadd $pos0 $d]

	# Append the new position.
	lappend traj $pos
	set pos0 $pos
    }

    return $traj
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

compute $name $regionName $trajTime $structPrefix $dcdFreq $dcdList $outDir $stride
exit
