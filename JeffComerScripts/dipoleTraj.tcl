# jcomer2@uiuc.edu

proc compute {name moiety structPrefix dcdList dcdFreq outDir startFrame} {
    set selText "all"
    set displayPeriod 1000
    set timestep 1.0
    set stride 1
    if {$stride <= 0} {set stride 1}
    set outFile $outDir/dipole_${name}.dat

    # Input:
    set psf $structPrefix.psf
    set pdb $structPrefix.pdb
    
    # Get the time change between frames in nanoseconds.
    set dt [expr 1.0e-6*$timestep*$dcdFreq*$stride]

    # Load the system.
    mol load psf $psf pdb $pdb
    set sel [atomselect top $selText]

    set out [open $outFile w]

    # Loop over the dcd files.
    set nFrames0 0
    foreach dcd $dcdList {
	# Load the trajectory.
	animate delete all
	mol addfile $dcd type dcd step $stride waitfor all
	set nFrames [molinfo top get numframes]
	puts [format "Reading %i frames." $nFrames]

	# Move forward computing at every step.
	for {set f $startFrame} {$f < $nFrames} {incr f} {
	    molinfo top set frame $f
	    set d [lindex [measure dipole $sel] 2]
	    set t [expr {$nFrames0 + $f*$dt}]

	    puts $out "$t $d"

	    # Write the time and distance.
	    if {$f % $displayPeriod == 0} {
		puts -nonewline [format "FRAME %i: " $f]
		puts "$t $d"
	    }
	}
	set nFrames0 [expr $nFrames+$nFrames0]
    }
    close $out

    $sel delete
    mol delete top
}
