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

    set z0 72.0
    set zMem 45.0
    set zBelow -65
    set expectedDisplace 6.0
    set poreRadius 25

    set xsc ../extra_layer_dopc_neg.restart.xsc
    set radIn [expr {$poreRadius - $expectedDisplace}]
    set radOut [expr {$poreRadius + $expectedDisplace}]

    set lowerSeg "LA65 LA66 LA69 LA70 LA71 LB37 LB39 LB40 LB41 LB42 LB46 LB47 LB48 LB51 LB52 LB53 LB57 LB58 LB59 LB63 LB64 LB65 LB66 LB69 LB70 LB71 LC37 LC39 LC41 LC42 LC46 LC47 LC48 LD51 LD52 LD53 LD55 LD56 LD57 LD58 LD59 LD60 LD61 LD62 LD63 LD64 LD65 LD66 LD67 LD68 LD69 LD70 LD71 LD72 LE37 LE38 LE39 LE40 LE41 LE42 LE43 LE44 LE45 LE46 LE47 LE48 LE49 LE50 LE51 LE52 LE53 LE54 LE55 LE56 LE57 LE58 LE59 LE60 LE61 LE62 LE63 LE64 LE65 LE66 LE67 LE68 LE69 LE70 LE71 LE72 LF37 LF38 LF39 LF40 LF41 LF42 LF43 LF44 LF45 LF46 LF47 LF48 LF49 LF50 LF51 LF52 LF53 LF55 LF56 LF57 LF58 LF64 LG61 LG62 LG67 LG68 LG69 LG72 LH38 LH39 LH43 LH44 LH45 LH46 LH49 LH50 LH51 LH52 LH54 LH55 LH56 LH57 LH60 LH61 LH62 LH63 LH67 LH68 LH69 LH72 LI38 LI39 LI43 LI44 LI45 LI49"

    # Define the regions.
    set rad0Sq [expr {$radIn*$radIn}]
    set rad1Sq [expr {$radOut*$radOut}]    
    set region(inTop) "resname PCGL and x^2+y^2<$rad0Sq and not segname $lowerSeg"
    set region(outTop) "resname PCGL and x^2+y^2>$rad1Sq and not segname $lowerSeg"
    set region(inBot) "resname PCGL and x^2+y^2<$rad0Sq and segname $lowerSeg"
    set region(outBot) "resname PCGL and x^2+y^2>$rad1Sq and segname $lowerSeg"
    set region(bot) "resname PCGL and segname $lowerSeg"
    set region(top) "resname PCGL and not segname $lowerSeg"
    set region(potFree) "name POT and z < $zBelow"
    set region(chlFree) "name CLA and z < $zBelow"
    set region(potLayer) "name POT and z < $z0 and z > $zMem and x^2+y^2>$rad1Sq"
    set region(chlLayer) "name CLA and z < $z0 and z > $zMem and x^2+y^2>$rad1Sq"
    set region(potLow) "name POT and z < -75"
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
