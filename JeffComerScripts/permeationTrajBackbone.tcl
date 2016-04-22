# jcomer2@uiuc.edu

proc compute {name moiety structPrefix dir dcdSet dcdFreq outDir} {
    set displayPeriod 200
    set stride 1
    set timestep 1.0
    set selText "name C1' H1' C2' H2' H2'' C3' O3' H3' C4' O4' H4' C5' O5' H5' H5'' O1P O2P P"

    set segA ADNA
    set segB BDNA

    # Input:
    set psf $structPrefix.psf
    set pdb $structPrefix.pdb
    set xsc $structPrefix.xsc
    set dcdPrefix $dir/nw_${name}
    set dcdSuffix .dcd
    # Output:
    set outPrefix $outDir/backbones_${name}.dat

    # Get the time change between frames in nanoseconds.
    set dt [expr 1.0e-6*$timestep*$dcdFreq*$stride]
    set dcdEnd [llength $dcdSet]
    
    # Read the system size from the xsc file.
    # Note: This only works for lattice vectors along the axes!
    set in [open $xsc r]
    foreach line [split [read $in] "\n"] {
	if {![string match "#*" $line]} {
	    set param [split $line]
	    puts $param
	    set lx [lindex $param 1]
	    set ly [lindex $param 5]
	    set lz [lindex $param 9]
	    break
	}
    }
    puts "NOTE: The system size is $lx $ly $lz.\n"
    close $in

    # Open the output file.
    set out [open $outPrefix w]

    # Load the system.
    mol load psf $psf pdb $pdb
    set selA [atomselect top "segname $segA"]
    set selB [atomselect top "segname $segB"]
    set resA [lsort -unique -integer [$selA get resid]]
    set resB [lsort -unique -integer [$selB get resid]]
    $selA delete
    $selB delete

    set resSelA {}
    set resSelB {}
    foreach a $resA b $resB {
	lappend resSelA [atomselect top "segname $segA and resid $a and ($selText)"]
	lappend resSelB [atomselect top "segname $segB and resid $b and ($selText)"]
    }

    # Loop over the dcd files.
    set nFrames0 0
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

	    set perm 0.0
	    foreach a $resSelA b $resSelB {
		set pos [$a get z]
		set n [$a num]
		set c 0
		foreach z $pos {
		    #if {$z > 0.5*$lz} {set z [expr $z - $lz]}
		    #if {$z <= -0.5*$lz} {set z [expr $z + $lz]}
		    if {$z < 0} {incr c}
		}
		set perm [expr $perm + 1.0*$c/$n]

		set pos [$b get z]
		set n [$b num]
		set c 0
		foreach z $pos {
		    #if {$z > 0.5*$lz} {set z [expr $z - $lz]}
		    #if {$z <= -0.5*$lz} {set z [expr $z + $lz]}
		    if {$z < 0} {incr c}
		}
		set perm [expr $perm + 1.0*$c/$n]
	    }

	    
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

    foreach a $resSelA b $resSelB {
	$a delete
	$b delete
    }
    mol delete top
}
