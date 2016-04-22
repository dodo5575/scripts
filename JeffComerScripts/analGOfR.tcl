# Author: jcomer2@illinois.edu

if {$argc < 8} {
    puts "$argv0 name srcMoiety destMoiety structPrefix outputDir dcdFreq stride startFrame gridFile dcdFile0 \[dcdFile1...\]"
    exit
}
set name [lindex $argv 0]
set srcMoiety [lindex $argv 1]
set destMoiety [lindex $argv 2]
set structPrefix [lindex $argv 3]
set outputDir [lindex $argv 4]
set dcdFreq [lindex $argv 5]
set stride [lindex $argv 6]
set startFrame [lindex $argv 7]
set gridFile [lindex $argv 8]
set dcdList [lrange $argv 9 end]

source $env(HOME)/scripts/vector.tcl
source $env(HOME)/scripts/gridForce.tcl

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

proc binGetVolSphere {binVar ir} {
    upvar $binVar bin

    set pi [expr {4.0*atan(1.0)}]
    set r0 [expr {$bin(x0) + $bin(dx)*$ir}]
    set r1 [expr {$bin(x0) + $bin(dx)*($ir+1)}]

    return [expr {4.0/3.0*$pi*($r1*$r1*$r1 - $r0*$r0*$r0)}]
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
proc binPut {binVar dataX} {
    upvar $binVar bin
 
    set ix [binGetCell bin $dataX]
    if {$ix >= 0 && $ix < $bin(n)} {
	incr bin(count,$ix)
	return 1
    }
    return 0
}
###############################
###############################

proc compute {name srcMoiety destMoiety structPrefix outputDir dcdFreq stride startFrame gridFile dcdList} {
    set displayPeriod 200
    set timestep 1.0
    set srcText "name $srcMoiety"
    set destText "name $destMoiety"
    if {$stride <= 0} {set stride 1}
    if {$startFrame < 1} {set startFrame 1}
    set bootNum 3

    # Load the grid for applying periodic boundaries.
    readDx grid $gridFile

    # Bin parameters.
    set r0 1.0
    set r1 19.0
    set nr 90

    # Make the main bin.
    binNew main $r0 $r1 $nr
    set bins [binGetSize main]

    # Make the boot bins.
    for {set i 0} {$i < $bootNum} {incr i} {
	binNew boot${i} $r0 $r1 $nr
	set bootSteps($i) 0
    }

    # Input:
    set pdb $structPrefix.pdb
    
    # Get the time change between frames in nanoseconds.
    set dt [expr {1.0e-6*$timestep*$dcdFreq*$stride}]

    # Load the system.
    mol load pdb $pdb
    set srcSel [atomselect top "$srcText"]
    set destSel [atomselect top "$destText"]
    set srcNum [$srcSel num]
    set destNum [$destSel num]
    
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
	    molinfo top set frame $f
	    # The index of the boot for this frame.
	    set bootInd [expr {int(floor(rand()*$bootNum))}]

	    # Get the positions at this step.
	    set srcList [$srcSel get {x y z}]
	    set destList [$destSel get {x y z}]

	    for {set i 0} {$i < $srcNum} {incr i} {
		for {set j 0} {$j < $destNum} {incr j} {
		    set r0 [lindex $srcList $i]
		    set r1 [lindex $destList $j]
		    set r [veclength [wrapDiff grid [vecsub $r1 $r0]]]
		    
		    # Add the data to the main bin set.
		    binPut main $r 
		    
		    # Add the data to the randomly chosen boot bin set.
		    binPut boot${bootInd} $r
		}
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

    $srcSel delete
    $destSel delete
    mol delete top

    if {$steps <= 1} {
	puts "NO DATA!"
	return
    }

    # The number of instances is equal to the number of steps times the number of src particles.
    # Get the main results.
    set quantList [computeQuant main final [expr {$steps*$srcNum}]]

    # Zero the error.
    foreach q $quantList {
	for {set is 0} {$is < $bins} {incr is} {
	    set err($is,$q) 0.0
	}
    }

    # Use the greatest deviation in the boot results to get the error.
    for {set b 0} {$b < $bootNum} {incr b} {
	computeQuant boot${b} bootQ [expr {$bootSteps($b)*$srcNum}]
	
	for {set is 0} {$is < $bins} {incr is} {
	    foreach q $quantList {
		set d [expr {abs($final($is,$q) - $bootQ($is,$q))}]
		if {$d > $err($is,$q)} { set err($is,$q) $d }
	    }
	}
    }
    
    # Write the results.
    set out [open $outputDir/gofr_${srcMoiety}-${destMoiety}_${name}.dat w]
    for {set is 0} {$is < $bins} {incr is} {
	# Write the position.
	puts -nonewline $out "$final($is,pos)"

	foreach q $quantList {
	    puts -nonewline $out " $final($is,$q) $err($is,$q)"
	}
	puts $out ""
    }
    close $out
   
    return
}

proc computeQuant {binVar quantVar steps} {
    upvar $binVar bin
    upvar $quantVar quant

    # Important constants.
    set concConvert 1660.5387; # particles/AA^3 -> mol/l

    set n [binGetSize bin]
    for {set ir 0} {$ir < $n} {incr ir} {
	# Get the bin position.
	set r [binGetCenter bin $ir]
	set vol [binGetVolSphere bin $ir]

	# Get the average number of particles.
	set count [binGetCount bin $ir]
	set nAvg [expr {double($count)/$steps}]

	# Compute the density in mol/L.
	set conc [expr {$concConvert*$nAvg/$vol}]
	
	# Return the important stuff.
	set quant($ir,pos) [expr {0.1*$r}]; # in nm
	set quant($ir,num) $nAvg; # in 1
	set quant($ir,conc) $conc; # in M
    }

    # Results are "(1)r (2)conc"

    return {conc}
}

compute $name $srcMoiety $destMoiety $structPrefix $outputDir $dcdFreq $stride $startFrame $gridFile $dcdList
exit
