# jcomer2@uiuc.edu

proc compute {name moiety structPrefix dcdList dcdFreq outDir startFrame} {
    set displayPeriod 200
    set timestep 1.0
    set stride 1
    set quantList {numK numCl}
    set textList {"name POT" "name CLA"}
    set regionTextList {"z >= -15 and  z < -5" "z >= -5 and z < 5" "z >= 5 and z < 15" "z >= 25 and z < 35"}
    set regionNameList {"z-10" "z0" "z10" "z30"}
    
    if {$stride <= 0} {set stride 1}
    if {$startFrame <= 0} {set startFrame 1}

    set s0 0.0
    set s1 20.0
    set ns 20
    set dz 10.0

    set avogadro 6.0221367e23
    
    set ds [expr {($s1-$s0)/$ns}]
    # Zero the sums.
    foreach quant $quantList {
	foreach reg $regionNameList {
	    puts "$reg $quant"

	    for {set is 0} {$is < $ns} {incr is} {
		set sum($quant,$reg,$is) 0.0
		set count($quant,$reg,$is) 0
	    }
	}
    }

    # Input:
    set psf $structPrefix.psf
    set pdb $structPrefix.pdb
  
    # Get the time change between frames in nanoseconds.
    set dt [expr 1.0e-6*$timestep*$dcdFreq*$stride]

    # Load the system.
    mol load pdb $pdb
    
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
	    # Get the time in nanoseconds for this frame.
	    #set t [expr {($nFrames0+$f)*$dt}]
	    incr totalCount

	    foreach reg $regionNameList regTex $regionTextList {
		foreach quant $quantList tex $textList {
		    # Get the current position.
		    molinfo top set frame $f
		    set sel [atomselect top "($tex) and ($regTex)"]
		    set num [$sel num]
		    set pos [$sel get {x y z}]
		    set ind [$sel get index]
		    $sel delete

		    if {$num < 1} { continue }
		    
		    # Get the last position.
		    molinfo top set frame [expr {$f-1}]
		    set sel [atomselect top "index $ind"]
		    set pos0 [$sel get {x y z}]
		    $sel delete

		    # Add each atom to the counter.
		    foreach r $pos r0 $pos0 {
			set m [vecscale 0.5 [vecadd $r $r0]]
			set d [vecsub $r $r0]
			set vz [expr {[lindex $d 2]/$dt}]; # in A/ns

			foreach {x y z} $m { break }
			set s [expr {sqrt($x*$x + $y*$y)}]
			set is [expr {int(floor(($s-$s0)/$ds))}]

			if {$is >= 0 && $is < $ns} {
			    set sum($quant,$reg,$is) [expr {$sum($quant,$reg,$is) + $vz}]
			    incr count($quant,$reg,$is)
			}
		    }
		}
	    }

	    # Write the time and distance.
	    if {$f % $displayPeriod == 0} {
		puts -nonewline [format "FRAME %i: " $f]
		puts "$totalCount"
	    }
	}
	set nFrames0 [expr $nFrames+$nFrames0]
    }

    # Write the results.
    set pi [expr {4.0*atan(1.0)}]
    foreach quant $quantList {
	foreach reg $regionNameList {
	    set out [open $outDir/${quant}_${name}_${reg}.dat w]
	    
	    for {set is 0} {$is < $ns} {incr is} {
		# Get the count.
		set n $count($quant,$reg,$is)
		if {$n <= 0} {
		    set vel 0.0
		} else {
		    set vel [expr {$sum($quant,$reg,$is)/double($n)}]
		}

		# Compute the velocity in A/ns.
		set s [expr {$s0 + $ds*$is}]
		set sn [expr {$s0 + $ds*($is+1)}]
		set sm [expr {0.5*($s+$sn)}]
		set rho $vel

		puts $out "$sm $rho"
	    }

	    close $out
	}
    }
    mol delete top
}
