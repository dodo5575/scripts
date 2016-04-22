 # jcomer2@uiuc.edu

proc compute {name moiety structPrefix dir dcdSet dcdFreq outDir} {
    set stride 1
    set timestep 1.0

    if {[string equal $moiety ecoRI]} {
	set selText0 "(segname ADNA and resid 557) or (segname BDNA and resid 612)"
	set selText1 "(segname ADNA and resid 567) or (segname BDNA and resid 602)"
    } elseif {[string equal $moiety bamHI]} {
	set selText0 "(segname ADNA and resid 34) or (segname BDNA and resid 29)"; 
	set selText1 "(segname ADNA and resid 44) or (segname BDNA and resid 19)";
    } else {
	puts "ERROR: Unrecognized moiety $moiety!"
    }

    # Input:
    set psf $structPrefix.psf
    set pdb $structPrefix.pdb
    set dcdPrefix $dir/nw_${name}
    set dcdSuffix .dcd
    # Output:
    set outPrefix $outDir/str_${name}.dat

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
    for {set dcd 0} {$dcd < $dcdEnd} {incr dcd} {
	# Load the trajectory.
	animate delete all
	set dcdNum [lindex $dcdSet $dcd]
	mol addfile "${dcdPrefix}${dcdNum}${dcdSuffix}" type dcd step $stride waitfor all
	set nFrames [molinfo top get numframes]
	puts [format "Reading %i frames." $nFrames]

	# Move forward computing the center-of-mass at every step.
	for {set f 0} {$f < $nFrames} {incr f} {
	    molinfo top set frame $f

	    # Get the time in nanoseconds for this frame.
	    set t [expr ($nFrames0+$f)*$dt]

	    set r0 [measure center $sel0 weight mass]
	    set r1 [measure center $sel1 weight mass]
	    # Get the distance in nanometers.
	    set d [expr 0.1*[veclength [vecsub $r1 $r0]]]

	    # Write the time and position.
	    puts $out "$t $d"
	    
	    puts -nonewline [format "FRAME %i: " $f]
	    puts "$t $d"
	}
	set nFrames0 [expr $nFrames+$nFrames0]
    }

    close $out
    $sel0 delete
    $sel1 delete
    mol delete top
}
