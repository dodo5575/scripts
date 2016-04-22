# jcomer2@uiuc.edu

proc compute {name moiety structPrefix dcdList dcdFreq outDir stride} {
    set displayPeriod 200
    set timestep 1.0
    set startFrame 0
    if {$stride < 1} {set stride 1}
    #set quantList {numK numCl numIon}
    #set textList {"name POT" "name CLA" "name POT CLA"}
    set quantList {numK}
    set textList {"name POT"}
    
    set regionTextList {}
    set regionNameList {}

    set z0 -35.0
    set z1 35.0
    set dz 0.5
    for {set z $z0} {$z < $z1} {set z [expr {$z + $dz}]} {
	set zn [expr {$z + $dz}]
	set zm [format "%g" [expr {0.5*($zn+$z)}]]
	lappend regionTextList "z >= $z and z < $zn"
	lappend regionNameList $zm
    }
    if {$stride <= 0} {set stride 1}

    # Init the sums.
    foreach quant $quantList {	
	foreach reg $regionNameList {
	    set sum($quant,$reg) 0.0
	    set sumSq($quant,$reg) 0.0
	}
    }

    # Input:
    set psf $structPrefix.psf
    set pdb $structPrefix.pdb
  
    # Get the time change between frames in nanoseconds.
    set dt [expr 1.0e-6*$timestep*$dcdFreq*$stride]

    # Load the system.
    mol load psf $psf pdb $pdb
    
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
	    molinfo top set frame $f

	    # Get the time in nanoseconds for this frame.
	    set t [expr {($nFrames0+$f)*$dt}]
	    incr totalCount

	    foreach quant $quantList tex $textList {
		foreach reg $regionNameList regTex $regionTextList {
		    set sel [atomselect top "($tex) and ($regTex)"]
		    set num [$sel num]
		    $sel delete
		    
		    set sum($quant,$reg) [expr {$sum($quant,$reg) + $num}]
		    set sumSq($quant,$reg) [expr {$sumSq($quant,$reg) + $num*$num}]
		}
	    }

	    # Write the time and distance.
	    if {$f % $displayPeriod == 0} {
		puts -nonewline [format "FRAME %i: " $f]
		puts "$t"
	    }
	}
	set nFrames0 [expr $nFrames+$nFrames0]
    }

    
    # Get the count.
    if {$totalCount == 0} { return }
    set n [expr {double($totalCount)}]

    # Write the results.
    foreach quant $quantList {
	set out [open $outDir/${quant}_${name}.dat w]
	
	foreach reg $regionNameList {
	    set mean [expr {$sum($quant,$reg)/$n}]
	    set err [expr {sqrt(($sum($quant,$reg) - $sumSq($quant,$reg)/$n)/($n-1)/$n)}]
	    puts $out "$reg $mean $err"
	}

	close $out
    }

    mol delete top
}
