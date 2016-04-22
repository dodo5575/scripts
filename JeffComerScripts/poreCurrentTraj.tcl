# jcomer2@uiuc.edu
 
# Return the current and number of charge carriers in the selection for the frame.
proc computeCurrent {frameCurr sel charge lz dt} {
    set frameLast [expr $frameCurr-1]
    
    # Get the current position of the selection.
    molinfo top set frame $frameCurr
    set z1 [$sel get z]
    
    # Get the last position of the selection.
    molinfo top set frame $frameLast
    set z0 [$sel get z]
    
    # Get the number of charge carriers.
    set num [$sel num]
    
    # Find the displacements in the z-direction and compute the current.
    set currentZ 0.0
    foreach a0 $z0 a1 $z1 {
	# Compensate for jumps across the periodic cell.
	set dz [expr $a1-$a0]
	if {[expr $dz > 0.5*$lz]} {set dz [expr $dz-$lz]}
	if {[expr $dz <-0.5*$lz]} {set dz [expr $dz+$lz]}

	# Compute the current in nanoamperes.
	set currentZ [expr $currentZ + $charge*$dz/($lz*$dt)*1.60217733e-1]
    }
    
    return [list $num $currentZ]
}

proc compute {name moiety structPrefix dir dcd dcdFreq outDir} {
    set displayPeriod 200
    set stride 1
    set timestep 1.0
    set poreHalfZ 40
    set selText [list "name POT and abs(z)<$poreHalfZ" "name CLA and abs(z)<$poreHalfZ"]
    set charge [list 1.0 -1.0]
    set nameList [list "K" "Cl"]
    set outName p${poreHalfZ}_curr

    # Input:
    set psf $structPrefix.psf
    set pdb $structPrefix.pdb
    set xsc $structPrefix.xsc
        
    # Get the time change between frames in nanoseconds.
    set dt [expr 1.0e-6*$timestep*$dcdFreq*$stride]

    # Read the system size from the xsc file.
    # Note: This only works for lattice vectors along the axes!
    set in [open $xsc r]
    foreach line [split [read $in] "\n"] {
	if {![string match "#*" $line]} {
	    set param [split $line]
	    puts $param
	    set lx [lindex $param 1]
	    set ly [lindex $param 5]
	    set lz [lindex $param 9]
	    break
	}
    }
    puts "NOTE: The system size is $lx $ly $lz.\n"
    close $in

    set poreZ [expr 2.0*$poreHalfZ]
    if {$poreZ > $lz} {
	puts "WARNING: Setting poreZ to $lz."
	set poreZ $lz
    }

    # Open the output files.
    set out {}
    foreach n $nameList {
	lappend out [open "${outDir}/${outName}${n}_${name}.dat" w]
    }
    set outTotal [open "${outDir}/${outName}_${name}.dat" w]

    # Load the system.
    mol load psf $psf pdb $pdb

    # Loop over the dcd files.
    set nFrames0 0
    foreach dcdFile $dcd {
	# Load the trajectory.
	animate delete all
	mol addfile $dir/$dcdFile type dcd step $stride waitfor all
	set nFrames [molinfo top get numframes]
	puts [format "Reading %i frames." $nFrames]

	# Move forward, computing
	# current at each step.
	for {set f 1} {$f < $nFrames} {incr f} {
	    # Get the time in nanoseconds for this frame.
	    set t [expr ($nFrames0+$f+0.5)*$dt]

	    # Compute the current for each selection.		
	    set currentZTotal 0.0
	    foreach s $selText q $charge o $out {
		#Write the number of carriers and the current for this selection.
		set sel [atomselect top $s]
		set data [computeCurrent $f $sel $q $poreZ $dt]
		set currentZ [lindex $data 1]
		puts $o "$t $currentZ"
		
		$sel delete
		set currentZTotal [expr $currentZTotal + $currentZ]
	    }
	    
	    # Write the total current.
	    puts $outTotal "$t $currentZTotal"
	    
	    # Update the display.
	    if {$f % $displayPeriod == 0} {
		puts -nonewline [format "FRAME %i: " $f]
		puts "$t $currentZTotal"
	    }
	}
	set nFrames0 [expr $nFrames+$nFrames0]
    }

    foreach o $out {
	close $o
    }
    close $outTotal
    mol delete top
}
