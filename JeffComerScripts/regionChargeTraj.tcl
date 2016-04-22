# jcomer2@uiuc.edu

proc compute {name moiety structPrefix dcdList dcdFreq outDir stride} {
    set displayPeriod 200
    set timestep 1.0
    set startFrame 0
    set quantList {dna}
    set textList {"segname ADNA BDNA"}

    set z0 -80
    set z1 80
    set dz 10

    # Read the system size from the xsc file.
    # Note: This only works for lattice vectors along the axes!
    set in [open $structPrefix.xsc r]
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

    # Set the region selection texts.
    # Handle periodic jumps.
    set regionTextList {}
    set regionNameList {}
    for {set z $z0} {$z < $z1} {set z [expr {$z+$dz}]} {
	set ze [expr {$z+$dz}]
	set zm [expr {0.5*($z+$ze)}]
	set zA [expr {$z + $lz}]
	set zeA [expr {$ze + $lz}]
	set zB [expr {$z - $lz}]
	set zeB [expr {$ze - $lz}]
	lappend regionTextList "(z >= $z and z < $ze) or (z >= $zA and z < $zeA) or (z >= $zB and z < $zeB)"
	lappend regionNameList "$zm"
    }

    # Open the output files.
    foreach quant $quantList {
	foreach reg $regionNameList {
	    set out($quant,$reg) [open $outDir/${quant}_${name}_${reg}.dat w]
	    puts "$reg $quant"
	}
    }

    # Input:
    set psf $structPrefix.psf
    set pdb $structPrefix.pdb
  
    # Get the time change between frames in nanoseconds.
    set dt [expr 1.0e-6*$timestep*$dcdFreq*$stride]

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

	    foreach reg $regionNameList regTex $regionTextList {
		foreach quant $quantList tex $textList {
		    set sel [atomselect top "($tex) and ($regTex)"]
		    #set num [$sel num]
		    set num [measure sumweights $sel weight charge]
		    $sel delete

		    puts $out($quant,$reg) "$t $num"
		}
	    }

	    # Write the time and distance.
	    if {$f % $displayPeriod == 0} {
		puts -nonewline [format "FRAME %i: " $f]
		puts "$t"
	    }
	}
	set nFrames0 [expr $nFrames+$nFrames0]
    }

    foreach quant $quantList {
	foreach reg $regionNameList {
	    close $out($quant,$reg)
	}
    }
    mol delete top
}
