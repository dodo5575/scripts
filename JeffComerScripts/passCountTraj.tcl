# jcomer2@uiuc.edu

proc compute {name moiety structPrefix dcdList dcdFreq outDir stride} {
    set planeZ 0
    set planeBuf 30
    set displayPeriod 200
    set timestep 1.0
    set selText "segname ADNA BDNA"
    set pre "passCount"
    if {$stride <= 0} {set stride 1}

    # Input:
    set psf $structPrefix.psf
    set pdb $structPrefix.pdb
    set xsc $structPrefix.xsc
    
    # Get the time change between frames in nanoseconds.
    set dt [expr 1.0e-6*$timestep*$dcdFreq*$stride]

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
    set out [open "${outDir}/${pre}_${name}.dat" w]

    # Load the system.
    mol load psf $psf
    # Load the first frame of the dcd.
    mol addfile [lindex $dcdList 0] type dcd first 0 last 0 waitfor all

    # Get each residue.
    set sel [atomselect top $selText]
    set resList [lsort -unique -integer [$sel get residue]]
    $sel delete

    set selList {}
    foreach res $resList {
	set s [atomselect top "residue $res"]
	lappend selList $s
	
	# Get the initial positions of each residue.
	set z [lindex [measure center $s weight mass] 2]
	set d [wrapDiffReal [expr {$z - $planeZ}] $lz]
	if {$d >= 0 && $d < $planeZ + $planeBuf} {
	    # in the top
	    set inTop($res) 1
	    set inBot($res) 0
	} elseif {$d < 0 && $d > $planeZ - $planeBuf} {
	    # in the bottom
	    set inTop($res) 0
	    set inBot($res) 1
	} else {
	    # elsewhere
	    set inTop($res) 0
	    set inBot($res) 0
	}
    }

    # The number passing through starts at zero.
    set passCount 0

    # Loop over the dcd files.
    set nFrames0 0
    foreach dcdFile $dcdList {
	# Load the trajectory.
	animate delete all
	mol addfile $dcdFile type dcd step $stride waitfor all
	set nFrames [molinfo top get numframes]
	puts [format "Reading %i frames." $nFrames]

	# Move forward, computing
	# current at each step.
	for {set f 0} {$f < $nFrames} {incr f} {
	    # Get the time in nanoseconds for this frame.
	    molinfo top set frame $f
	    set t [expr ($nFrames0+$f+0.5)*$dt]

	    # Compute the current for each selection.
	    foreach s $selList res $resList  {
		set z [lindex [measure center $s weight mass] 2]
		set d [wrapDiffReal [expr {$z - $planeZ}] $lz]
		if {$d >= 0 && $d < $planeZ + $planeBuf} {
		    # In the top.
		    # Did it pass through the plane in the negative direction?
		    if {$inBot($res)} {
			#puts pass-
			incr passCount -1
		    }

		    set inTop($res) 1
		    set inBot($res) 0
		} elseif {$d < 0 && $d > $planeZ - $planeBuf} {
		    # In the bottom.
		    # Did it pass through the plane in the positive direction?
		    if {$inTop($res)} {
			#puts pass+
			incr passCount 1
		    }

		    set inTop($res) 0
		    set inBot($res) 1
		} else {
		    # Elsewhere
		    set inTop($res) 0
		    set inBot($res) 0
		}
	    }
	    puts $out "$t $passCount"

	    # Update the display.
	    if {$f % $displayPeriod == 0} {
		puts -nonewline [format "FRAME %i: " $f]
		puts "$t $passCount"
	    }
	}
	set nFrames0 [expr $nFrames+$nFrames0]
    }


    foreach s $selList {
	$s delete
    }
    mol delete top
}

proc wrapDiffReal {x l} {
    set l [expr {double($l)}]
    set image [expr {int(floor($x/$l))}]
    set x [expr {$x - $image*$l}]

    if {$x >= 0.5*$l} { set x [expr {$x - $l}] }
    return $x
}
