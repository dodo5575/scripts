# Author: jcomer2@illinois.edu

if {$argc < 8} {
    puts "$argv0 name moiety structPrefix outputDir dcdFreq stride startFrame regionName dcdFile0 \[dcdFile1...\]"
    exit
}
set name [lindex $argv 0]
set moiety [lindex $argv 1]
set structPrefix [lindex $argv 2]
set outputDir [lindex $argv 3]
set dcdFreq [lindex $argv 4]
set stride [lindex $argv 5]
set startFrame [lindex $argv 6]
set regionName [lindex $argv 7]
set dcdList [lrange $argv 8 end]

###############################
###############################
# A bin array-object
proc binNew {binVar x0 x1 n} {
    upvar $binVar bin

    set bin(dx) [expr {double($x1-$x0)/$n}]
    set bin(x0) $x0
    set bin(n) $n

    # Zero the cells.
    for {set ix 0} {$ix < $bin(n)} {incr ix} {
	set bin(count,$ix) 0
	set bin(y,$ix) 0.0
	set bin(y2,$ix) 0.0
    }
    return $bin(dx)
}

proc binGetCell {binVar x} {
    upvar $binVar bin
    set ix [expr {int(floor(($x-$bin(x0))/$bin(dx)))}]
    return $ix
}

proc binGetCellCyl {binVar pos} {
    upvar $binVar bin

    foreach {x y z} $pos { break }
    set s [expr {sqrt($x*$x + $y*$y)}]
    return [binGetCell bin $s]
}

proc binGetCenter {binVar ix} {
    upvar $binVar bin

    return [expr {$bin(x0) + $bin(dx)*($ix + 0.5)}]
}

proc binGetAreaCyl {binVar is} {
    upvar $binVar bin

    set pi [expr {4.0*atan(1.0)}]
    set s0 [expr {$bin(x0) + $bin(dx)*$is}]
    set s1 [expr {$bin(x0) + $bin(dx)*($is+1)}]

    return [expr {$pi*($s1*$s1 - $s0*$s0)}]
}

proc binGetCount {binVar ix} {
    upvar $binVar bin
    
    return $bin(count,$ix)
}

proc binGetValue {binVar ix} {
    upvar $binVar bin
    
    return $bin(y,$ix)
}

proc binGetValue2 {binVar ix} {
    upvar $binVar bin
    
    return $bin(y2,$ix)
}

proc binGetSize {binVar} {
    upvar $binVar bin

    return $bin(n)
}

# Add some data to the bins.
proc binPut {binVar dataX dataY} {
    upvar $binVar bin
 
    set ix [binGetCell bin $dataX]
    if {$ix >= 0 && $ix < $bin(n)} {
	incr bin(count,$ix)
	set bin(y,$ix) [expr {$bin(y,$ix) + $dataY}]
	set bin(y2,$ix) [expr {$bin(y2,$ix) + $dataY*$dataY}]
	return 1
    }
    return 0
}
###############################
###############################

proc compute {name moiety structPrefix outputDir dcdFreq stride startFrame regionName dcdList} {
    set displayPeriod 1000
    set timestep 1.0
    set selText "name $moiety"
    if {$stride <= 0} {set stride 1}
    if {$startFrame < 1} {set startFrame 1}
    set bootNum 3

    # Bin parameters.
    set s0 0.0
    set s1 20.0
    set ns 60

    # Get the region text.
    set regText(mid) {"abs(z) < 12" 24.0}
    ##
    set regionText [lindex $regText($regionName) 0]
    set regionLength [lindex $regText($regionName) 1]

    # Set the ion's charge
    set ionCharge(SOD) 1.0.
    set ionCharge(POT) 1.0
    set ionCharge(CLA) -1.0
    ##
    set charge $ionCharge($moiety)

    # Make the main bin.
    binNew main $s0 $s1 $ns
    set bins [binGetSize main]

    # Make the boot bins.
    for {set i 0} {$i < $bootNum} {incr i} {
	binNew boot${i} $s0 $s1 $ns
	set bootSteps($i) 0
    }

    # Input:
    set pdb $structPrefix.pdb
    
    # Get the time change between frames in nanoseconds.
    set dt [expr {1.0e-6*$timestep*$dcdFreq*$stride}]

    # Load the system.
    mol load pdb $pdb
    
    # Loop over the dcd files.
    set steps 0
    foreach dcd $dcdList {
	# Load the trajectory.
	animate delete all
	mol addfile $dcd type dcd step $stride waitfor all
	set nFrames [molinfo top get numframes]
	puts [format "Reading %i frames." $nFrames]

	# Move forward computing at every step.
	for {set f $startFrame} {$f < $nFrames} {incr f} {
	    set f0 [expr {$f-1}]

	    # The index of the boot for this frame.
	    set bootInd [expr {int(floor(rand()*$bootNum))}]

	    # Get the positions at this step and last.
	    set sel [atomselect top "($selText) and ($regionText)" frame $f]
	    set pos1 [$sel get {x y z}]
	    $sel frame $f0
	    set pos0 [$sel get {x y z}]
	    $sel delete

	    # Add the data from each ion.
	    foreach r0 $pos0 r1 $pos1 {
		foreach {x y z} $r0 { break }
		set s [expr {sqrt($x*$x + $y*$y)}]
		set dz [expr {[lindex $r1 2]-[lindex $r0 2]}]

		# Add the data to the main bin set.
		binPut main $s $dz
	
		# Add the data to the randomly chosen boot bin set.
		binPut boot${bootInd} $s $dz
	    }

	    # Write the status.
	    if {$f % $displayPeriod == 0} {
		puts -nonewline [format "FRAME %i: " $f]
		puts "$steps"
	    }

	    # Count the valid frames in each bin set.
	    incr steps
	    incr bootSteps($bootInd)
	}; # End frame loop.
    }

    if {$steps <= 1} {
	puts "NO DATA!"
	return
    }
    
    # Get the main results.
    set quantList [computeQuant main final $steps $dt $regionLength $charge]

    # Zero the error.
    foreach q $quantList {
	for {set is 0} {$is < $bins} {incr is} {
	    set err($is,$q) 0.0
	}
    }

    # Use the greatest deviation in the boot results to get the error.
    for {set b 0} {$b < $bootNum} {incr b} {
	computeQuant boot${b} bootQ $bootSteps($b) $dt $regionLength $charge
	
	for {set is 0} {$is < $bins} {incr is} {
	    foreach q $quantList {
		set d [expr {abs($final($is,$q) - $bootQ($is,$q))}]
		if {$d > $err($is,$q)} { set err($is,$q) $d }
	    }
	}
    }
    
    # Write the results.
    set out [open $outputDir/cyl_${moiety}_${name}.dat w]
    for {set is 0} {$is < $bins} {incr is} {
	# Write the position.
	puts -nonewline $out "$final($is,pos)"

	foreach q $quantList {
	    puts -nonewline $out " $final($is,$q) $err($is,$q)"
	}
	puts $out ""
    }
    close $out

    mol delete top
    return
}

proc computeQuant {binVar quantVar steps dt regionLength charge} {
    upvar $binVar bin
    upvar $quantVar quant

    # Important constants.
    set concConvert 1660.5387; # particles/AA^3 -> mol/l
    set currConvert 160.21765; # e/ns -> pA

    set n [binGetSize bin]
    for {set is 0} {$is < $n} {incr is} {
	# Get the bin position.
	set s [binGetCenter bin $is]
	set area [binGetAreaCyl bin $is]
	set vol [expr {$area*$regionLength}]

	# Get the average number of particles.
	set count [binGetCount bin $is]
	set nAvg [expr {double($count)/$steps}]

	# Get the average displacements.
	if {$count > 0} {
	    set dispAvg [expr {[binGetValue bin $is]/$steps}]
	    set disp2Avg [expr {[binGetValue2 bin $is]/$steps}]
	    set dPartAvg [expr {$dispAvg/$nAvg}]
	    set d2PartAvg [expr {$disp2Avg/$nAvg}]
	} else {
	    set dispAvg 0.0
	    set disp2Avg 0.0
	    set dPartAvg 0.0
	    set d2PartAvg 0.0
	}

	# Get the variance of the position.
	set var [expr {$d2PartAvg - $dPartAvg*$dPartAvg}]
	set diffuse [expr {$var/(2.0*$dt)}]

	# Compute the density in mol/L.
	set conc [expr {$concConvert*$nAvg/$vol}]

	# Compute the velocity in AA/ns and the current in pA.
	set vel [expr {$dPartAvg/$dt}]
	set pre [expr {$currConvert*$charge/($dt*$regionLength)}]
	set curr [expr {$pre*$dispAvg}]
	set currDensity [expr {$curr/$area}]
	
	# Return the important stuff.
	set quant($is,pos) [expr {0.1*$s}]; # in nm
	set quant($is,num) $nAvg; # in 1
	set quant($is,disp) [expr {0.1*$dispAvg}]; # in nm
	set quant($is,vel) [expr {0.1*$vel}]; # in nm/ns
	set quant($is,conc) $conc; # in M
	set quant($is,diffuse) [expr {0.01*$diffuse}]; # in nm^2/ns
	set quant($is,curr) $curr; # in pA
	set quant($is,j) [expr {$currDensity*100.0}]; # in pA/nm^2
    }

    # Results are "(1)s (2)n (3)en (4)d (5)ed (6)v (7)ev (8)c (9)ec (10)D (11)eD (12)I (13)eI (14)j (15)ej"

    return {num disp vel conc diffuse curr j}
}

compute $name $moiety $structPrefix $outputDir $dcdFreq $stride $startFrame $regionName $dcdList
exit
