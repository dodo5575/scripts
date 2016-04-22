# jcomer2@uiuc.edu

proc compute {name moiety structPrefix dcdSet dcdFreq outDir startFrame} {
    set displayPeriod 200
    set stride 1
    set timestep 1.0

    if {[string equal $moiety ecoRI]} {
	set selText0 "segname BDNA and resid 620 and name P"; # cognate sequence
	set selText1 "segname PRT2 and resid 113"; # inner recognition helices
    } elseif {[string equal $moiety bamHI]} {
	set selText0 "segname BDNA and resid 36 and name P"; # reactive P
	set selText1 "segname PRTA and resid 113"; # E113
    } else {
	puts "ERROR: Unrecognized moiety $moiety!"
    }

    # Input:
    set psf $structPrefix.psf
    set pdb $structPrefix.pdb
    set dcdSuffix .dcd
    # Output:
    set outPrefix $outDir/dist_${name}.dat

    # Get the time change between frames in nanoseconds.
    set dt [expr 1.0e-6*$timestep*$dcdFreq*$stride]
    set dcdEnd [llength $dcdSet]

    # Open the output file.
    set out [open $outPrefix w]

    # Load the system.
    mol load psf $psf pdb $pdb
    set sel0 [atomselect top $selText0]
    set sel1 [atomselect top $selText1]
    puts "[$sel0 num] atoms and [$sel1 num] atoms."

    # Loop over the dcd files.
    set nFrames0 1
    foreach dcd $dcdSet {
	# Load the trajectory.
	animate delete all
	mol addfile $dcd type dcd step $stride waitfor all
	set nFrames [molinfo top get numframes]
	puts [format "Reading %i frames." $nFrames]

	# Move forward computing the center-of-mass at every step.
	for {set f $startFrame} {$f < $nFrames} {incr f} {
	    molinfo top set frame $f

	    # Get the time in nanoseconds for this frame.
	    set t [expr ($nFrames0+$f)*$dt]

	    set r0 [measure center $sel0 weight mass]
	    set r1 [measure center $sel1 weight mass]
	    # Get the distance in nanometers.
	    set d [expr 0.1*[veclength [vecsub $r1 $r0]]]
	    
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
    $sel0 delete
    $sel1 delete
    mol delete top
}
