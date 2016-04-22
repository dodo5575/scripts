# jcomer2@uiuc.edu

proc compute {name moiety structPrefix dcdList dcdFreq outDir startFrame} {
    set displayPeriod 200
    set stride 1
    set timestep 1.0
    set centerResidueRadius 5
    set selText "segname ADNA BDNA"
    # Input:
    set psf $structPrefix.psf
    set pdb $structPrefix.pdb
    # Output:
    set outPrefix $outDir/bp_${name}.dat

    set cutoff 4
    set bondA "resname ADE and name N1"
    set bondT "resname THY and name H3"
    set bondC "resname CYT and name N3"
    set bondG "resname GUA and name H1"
    set pairAT "($selText and $bondA) and within $cutoff of ($bondT)"
    set pairCG "($selText and $bondC) and within $cutoff of ($bondG)"

    # Get the time change between frames in nanoseconds.
    set dt [expr 1.0e-6*$timestep*$dcdFreq*$stride]

    # Open the output file.
    set out [open $outPrefix w]

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

	    set numPairs 0
	    set selAT [atomselect top $pairAT]
	    set selCG [atomselect top $pairCG]
	    
	    set numPairs [expr $numPairs + [$selAT num]]
	    set numPairs [expr $numPairs + [$selCG num]]
	    $selAT delete
	    $selCG delete

	    # Write the time and distance.
	    puts $out "$t $numPairs"
	    if {$f % $displayPeriod == 0} {
		puts -nonewline [format "FRAME %i: " $f]
		puts "$t $numPairs"
	    }
	}
	set nFrames0 [expr $nFrames+$nFrames0]
    }

    close $out
    
    foreach s $selList {
	$s delete
    }

    mol delete top
}
