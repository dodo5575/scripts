# jcomer2@uiuc.edu

proc compute {name moiety structPrefix dcdList dcdFreq outDir stride} {
    set displayPeriod 200
    set timestep 1.0
    set startFrame 0
    set z1 53
    set z0 23
    set vol [expr {100.0*100.0*($z1 - $z0)}]
    set selText "name OH2"
    set unifiedMass 1.660538782e-27

    # Input:
    set psf $structPrefix.psf
    set pdb $structPrefix.pdb
    # Output:
    set outPrefix $outDir/density_${name}.dat

    # Get the time change between frames in nanoseconds.
    set dt [expr 1.0e-6*$timestep*$dcdFreq*$stride]

    # Open the output file.
    set out [open $outPrefix w]

    # Load the system.
    mol load psf $psf pdb $pdb
    
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
	    set sel [atomselect top "($selText) and z >= $z0 and z < $z1"]
	    #set d [$sel num]
	    set m [measure sumweights $sel weight mass]
	    $sel delete
	    # Get the density in kg/m^3.
	    set d [expr {$unifiedMass*1.0e30*$m/$vol}]

	    # Write the time and distance.
	    puts $out "$t $d"
	    if {$f % $displayPeriod == 0} {
		puts -nonewline [format "FRAME %i: " $f]
		puts "$t $d"
	    }
	}
	set nFrames0 [expr $nFrames+$nFrames0]
    }

    close $out
    mol delete top
}
