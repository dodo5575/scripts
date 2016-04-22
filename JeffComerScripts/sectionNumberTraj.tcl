# jcomer2@uiuc.edu

proc compute {name moiety structPrefix dcdList dcdFreq outDir startFrame} {
    set displayPeriod 1000
    set timestep 1.0
    set stride 1
    set quantList {numK numCl numIon}
    set textList {"name POT" "name CLA" "name POT CLA"}
    if {$stride <= 0} {set stride 1}

    set z0 -35.0
    set z1 35.0
    set nz 70
    
    set dz [expr {($z1-$z0)/$nz}]
    # Zero the sums.
    foreach quant $quantList {
	puts "$quant"

	for {set iz 0} {$iz < $nz} {incr iz} {
	    set sum($quant,$iz) 0
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
    set count 0
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
	    #set t [expr {($nFrames0+$f)*$dt}]
	    incr count

	    foreach quant $quantList tex $textList {
		set sel [atomselect top "($tex)"]
		set pos [$sel get {x y z}]
		$sel delete

		# Add each atom to the counter.
		foreach r $pos {
		    foreach {x y z} $r { break }
		    set iz [expr {int(floor(($z-$z0)/$dz))}]

		    if {$iz >= 0 && $iz < $nz} {
			incr sum($quant,$iz)
		    }
		}
	    }

	    # Write the time and distance.
	    if {$f % $displayPeriod == 0} {
		puts -nonewline [format "FRAME %i: " $f]
		puts "$count"
	    }
	}
	set nFrames0 [expr $nFrames+$nFrames0]
    }

    # Write the results.
    set pi [expr {4.0*atan(1.0)}]
    foreach quant $quantList {
	set out [open $outDir/${quant}_${name}.dat w]
	
	for {set iz 0} {$iz < $nz} {incr iz} {
	    # Get the count.
	    set n [expr {double($sum($quant,$iz))/$count}]

	    # Compute the density in mol/L.
	    set z [expr {$z0 + $dz*$iz}]
	    set zn [expr {$z0 + $dz*($iz+1)}]
	    set zm [expr {0.5*($z+$zn)}]
	    puts $out "$zm $n"
	}
	close $out
    }
    mol delete top
}
