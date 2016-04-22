# jcomer2@uiuc.edu

proc compute {name moiety structPrefix dir dcdSet dcdFreq outDir} {
    set stride 1
    set timestep 1.0
    set baseStep 4
    set displayPeriod 200

    if {[string equal $moiety ecoRI]} {
	set baseResA 523
	set baseResB 646
	set baseN 62
	set segNameA ADNA
	set segNameB BDNA
    } elseif {[string equal $moiety bamHI]} {
	set baseResA 1
	set baseResB 62
	set baseN 62
	set segNameA ADNA
	set segNameB BDNA
    } else {
	puts "ERROR: Unrecognized moiety $moiety!"
    }

    # Input:
    set psf $structPrefix.psf
    set pdb $structPrefix.pdb
    set dcdPrefix $dir/nw_${name}
    set dcdSuffix .dcd
    # Output:
    set outPrefix $outDir/strain_${name}.dat

    # Get the time change between frames in nanoseconds.
    set dt [expr 1.0e-6*$timestep*$dcdFreq*$stride]
    set dcdEnd [llength $dcdSet]

    # Open the output file.
    set out [open $outPrefix w]

    # Load the system.
    mol load psf $psf pdb $pdb

    set selList0 {}
    set selList1 {}
    set last 0
    for {set i $baseStep} {$i < $baseN} {incr i $baseStep} {
	set resA0 [expr $baseResA + $last]
	set resB0 [expr $baseResB - $last]
	set resA [expr $baseResA + $i]
	set resB [expr $baseResB - $i]
	
	set s0 [atomselect top "(segname $segNameA and resid $resA0) or (segname $segNameB and resid $resB0)"]
	set s1 [atomselect top "(segname $segNameA and resid $resA) or (segname $segNameB and resid $resB)"]
	
	if {[$s0 num] > 0 && [$s1 num]} {
	    lappend selList0 $s0
	    lappend selList1 $s1
	} else {
	    puts "ERROR: No DNA for selection ${segNameA}:${resA} ${segNameB}:${resB}" 
	    exit
	}

	set last $i
    }

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

	    foreach s0 $selList0 s1 $selList1 {
		set r0 [measure center $s0 weight mass]
		set r1 [measure center $s1 weight mass]
		set z [expr 0.1*0.5*([lindex $r0 2] + [lindex $r1 2])]
		# Get the distance in nanometers per basepair.
		set d [expr 0.1*[veclength [vecsub $r1 $r0]]/$baseStep]

		# Write the time and position and strain.
		puts $out "$t $z $d"
	    }
	    
	    if {$f % $displayPeriod == 0} {
		puts -nonewline [format "FRAME %i: " $f]
		puts "$t $d"
	    }
	}
	set nFrames0 [expr $nFrames+$nFrames0]
    }

    close $out
    foreach s0 $selList0 s1 $selList1 {
	$s0 delete
	$s1 delete
    }
    mol delete top
}
