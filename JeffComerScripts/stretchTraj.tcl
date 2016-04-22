# jcomer2@uiuc.edu

proc compute {name moiety structPrefix dcdList dcdFreq outDir startFrame} {
    set displayPeriod 200
    set stride 1
    set timestep 1.0
    set selText0 "segname ADNA and resid 1"
    set selText1 "segname ADNA and resid 20"
    set num 19
    set backbone "C1' H1' C2' H2' H2'' C3' O3' H3' C4' O4' H4' C5' O5' H5' H5'' O1P O2P P"

    # Input:
    set psf $structPrefix.psf
    set pdb $structPrefix.pdb
    # Output:
    set outPrefix $outDir/str_${name}.dat

    # Get the time change between frames in nanoseconds.
    set dt [expr 1.0e-6*$timestep*$dcdFreq*$stride]

    # Open the output file.
    set out [open $outPrefix w]

    # Load the system.
    mol load psf $psf pdb $pdb
    set sel0 [atomselect top $selText0]
    set sel1 [atomselect top $selText1]
    
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

	    set r0 [measure center $sel0 weight mass]
	    set r1 [measure center $sel1 weight mass]
	    set d [expr {[veclength [vecsub $r1 $r0]]/$num}]
	    
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
