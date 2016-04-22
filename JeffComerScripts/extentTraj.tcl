# jcomer2@illinois.edu

proc compute {name moiety structPrefix dcdList dcdFreq outDir startFrame} {
    set displayPeriod 200
    set stride 1
    set timestep 1.0
    set selText "segname ADNA BDNA"

    # Input:
    set psf $structPrefix.psf
    set pdb $structPrefix.pdb
    # Output:
    set outPrefix $outDir/extent_${name}.dat

    # Get the time change between frames in nanoseconds.
    set dt [expr 1.0e-6*$timestep*$dcdFreq*$stride]

    # Open the output file.
    set out [open $outPrefix w]

    # Load the system.
    mol load psf $psf pdb $pdb
    set sel [atomselect top $selText]
    
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

	    # Get the time in nanoseconds for this frame.
	    set t [expr {($nFrames0+$f)*$dt}]

	    #set r [lindex [measure center $sel weight mass] 2]
	    set zList [$sel get z]
	    set minZ [lindex $zList 0]
	    foreach z $zList {
		if {$z < $minZ} { set minZ $z }
	    }
	    set r [expr {0.1*$minZ}]

	    # Write the time and distance.
	    puts $out "$t $r"
	    if {$f % $displayPeriod == 0} {
		puts -nonewline [format "FRAME %i: " $f]
		puts "$t $r"
	    }
	}
	set nFrames0 [expr $nFrames+$nFrames0]
    }

    close $out

    $sel delete
    mol delete top
}
