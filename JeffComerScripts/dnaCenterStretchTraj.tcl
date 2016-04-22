# jcomer2@uiuc.edu

proc compute {name moiety structPrefix dcdList dcdFreq outDir startFrame} {
    set displayPeriod 200
    set stride 1
    set timestep 1.0
    set centerResidueRadius 5
    set selText "segname ADNA"
    # Input:
    set psf $structPrefix.psf
    set pdb $structPrefix.pdb
    # Output:
    set outPrefix $outDir/len_${name}.dat
    set outPrefixInverse $outDir/invlen_${name}.dat

    # Get the time change between frames in nanoseconds.
    set dt [expr 1.0e-6*$timestep*$dcdFreq*$stride]

    # Open the output file.
    set out [open $outPrefix w]
    set outInverse [open $outPrefixInverse w]

    # Load the system.
    mol load psf $psf pdb $pdb
    set sel [atomselect top $selText]
    set resList [lsort -unique -integer [$sel get resid]]
    $sel delete

    # Select the residues.
    set selList {}
    foreach res $resList {
	lappend selList [atomselect top "($selText) and resid $res"]
    }

    # Loop over the dcd files.
    set nFrames0 0
    foreach dcd $dcdList {
	# Load the trajectory.
	animate delete all
	mol addfile $dcd type dcd step $stride waitfor all
	set nFrames [molinfo top get numframes]
	puts [format "Reading %i frames." $nFrames]

	# Move forward computing at every step.
	for {set f 0} {$f < $nFrames} {incr f} {
	    molinfo top set frame $f

	    # Get the time in nanoseconds for this frame.
	    set t [expr ($nFrames0+$f)*$dt]

	    # Get the residue positions and find the center residue.
	    set resCen [lindex $resList 0]
	    set resCenZ 1e20
	    foreach s $selList res $resList {
		set r [measure center $s weight mass]
		set z [lindex $r 2]

		if {abs($z) < $resCenZ} {
		    set resCenZ [expr abs($z)]
		    set resCen $res
		}
		set pos($res) $z
	    }
	    
	    # Make sure the residue range is valid.
	    set res0 [expr $resCen-$centerResidueRadius]
	    set res1 [expr $resCen+$centerResidueRadius]
	    if {$res0 < [lindex $resList 0]} {
		set res0 [lindex $resList 0]
	    }
	    if {$res1 < [lindex $resList 0]} {
		set res1 [lindex $resList 0]
	    }
	    if {$res0 > [lindex $resList end]} {
		set res0 [lindex $resList end]
	    }
	    if {$res1 > [lindex $resList end]} {
		set res1 [lindex $resList end]
	    }

	    if {$res0 != $res1} {
		set l [expr abs(0.1*($pos($res1)-$pos($res0))/($res1-$res0))]
		set il [expr 1.0/$l]
	    } else {
		set l 0.0
	    }
	    #puts "$pos($res1) $pos($res0)"

	    # Write the time and distance.
	    puts $out "$t $l"
	    puts $outInverse "$t $il"
	    if {$f % $displayPeriod == 0} {
		puts -nonewline [format "FRAME %i: " $f]
		puts "$t $l $il"
	    }
	}
	set nFrames0 [expr $nFrames+$nFrames0]
    }

    close $out
    close $outInverse
    
    foreach s $selList {
	$s delete
    }

    mol delete top
}
