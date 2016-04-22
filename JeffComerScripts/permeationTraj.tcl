# jcomer2@uiuc.edu

proc compute {name moiety structPrefix dcdList dcdFreq outDir startFrame} {
    set stride 1
    set displayPeriod 200
    set timestep 2.0
    set selText "segname ADNA BDNA"
    # Input:
    set psf $structPrefix.psf
    set pdb $structPrefix.pdb
    # Output:
    set outPrefix $outDir/nucleotides_${name}.dat

    # Get the time change between frames in nanoseconds.
    set dt [expr 1.0e-6*$timestep*$dcdFreq*$stride]

    # Open the output file.
    set out [open $outPrefix w]

    # Load the system.
    mol load psf $psf pdb $pdb
    set sel [atomselect top $selText]
    set res [lsort -unique [$sel get {segname resid}]]
    $sel delete

    set resSel {}
    foreach r $res {
	lappend resSel [atomselect top "segname [lindex $r 0] and resid [lindex $r 1]"]
    }

    # Loop over the dcd files.
    set nFrames0 0
    foreach dcd $dcdList {
	# Load the trajectory.
	animate delete all
	mol addfile $dcd type dcd step $stride waitfor all
	set nFrames [molinfo top get numframes]
	puts [format "Reading %i frames." $nFrames]

	if {$nFrames0 == 0} {
	    molinfo top set frame 0

	    set perm 0.0
	    foreach r $resSel {
		set c 0
		set pos [$r get z]
		set n [$r num]
		foreach z $pos {
		    if {$z < 0} {incr c}
		}
		set perm [expr {$perm + 1.0*$c/$n}]
	    }
	    
	    set perm0 $perm
	}

	# Move forward computing at every step.
	for {set f 0} {$f < $nFrames} {incr f} {
	    molinfo top set frame $f

	    # Get the time in nanoseconds for this frame.
	    set t [expr ($nFrames0+$f)*$dt]

	    # Compute the number of permeated atoms for each residue.
	    set perm 0.0
	    foreach r $resSel {
		set c 0
		set pos [$r get z]
		set n [$r num]
		foreach z $pos {
		    if {$z < 0} {incr c}
		}
		set perm [expr {$perm + 1.0*$c/$n}]
	    }
	    set perm [expr {$perm - $perm0}]
	    
	    # Write the time and distance.
	    puts $out "$t $perm"
	    if {$f % $displayPeriod == 0} {
		puts -nonewline [format "FRAME %i: " $f]
		puts "$t $perm"
	    }
	}
	set nFrames0 [expr $nFrames+$nFrames0]
    }

    close $out

    foreach r $resSel {
	$r delete
    }
    mol delete top
}
