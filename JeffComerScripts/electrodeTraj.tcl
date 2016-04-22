# Author: Jeff Comer <jcomer2@illinois.edu>
# Write the potential (mV) on an electrode approximated by probePdb.

proc compute {name moiety structPrefix dir dcdSet dcdFreq outDir probePdb probeName} {
    set displayPeriod 200
    set stride 1
    set timestep 1.0
    set selText [list "segname ADNA BDNA"]
    set conc 0.1
    # Input:
    set psf $structPrefix.psf
    set pdb $structPrefix.pdb
    set xsc $structPrefix.xsc
    set dcdPrefix $dir/nw_${name}
    set dcdSuffix .dcd
    
    # Get the time change between frames in nanoseconds.
    set dt [expr 1.0e-6*$timestep*$dcdFreq*$stride]
    set dcdEnd [llength $dcdSet]

    # Calcuations
    set elemCharge 1.602176487e-19
    set eps0 8.854187817e-12
    set pi [expr 4.0*atan(1.0)]
    set debye [expr 0.3163*sqrt($conc)]
    set eps 80.0
    set elecConst [expr $elemCharge/(4.0*$pi*$eps0*$eps)]

    # Read the system size from the xsc file.
    # Note: This only works for lattice vectors along the axes!
    set in [open $xsc r]
    foreach line [split [read $in] "\n"] {
	if {![string match "#*" $line]} {
	    set param [split $line]
	    puts $param
	    set cellX [lrange $param 1 3]
	    set cellY [lrange $param 4 6]
	    set cellZ [lrange $param 7 9]
	    break
	}
    }
    puts "NOTE: The system size is $cellX, $cellY, $cellZ.\n"
    close $in

    # Open the output files.
    set out [open "${outDir}/${probeName}_${name}.dat" w]

    # Load the probes.
    set probeMol [mol new $probePdb]
    set probeSel [atomselect $probeMol all]
    set probeInd [$probeSel get index]
    set probeNum [llength $probeInd]
    $probeSel delete

    # Load the system.
    set mole [mol load psf $psf pdb $pdb]
    set sel [atomselect $mole $selText]
    set charge [$sel get charge]
    set indList [$sel get index]

    # Loop over the dcd files.
    set nFrames0 0
    for {set dcd 0} {$dcd < $dcdEnd} {incr dcd} {
	# Load the trajectory.
	animate delete all
	set dcdNum [lindex $dcdSet $dcd]
	mol addfile "${dcdPrefix}${dcdNum}${dcdSuffix}" type dcd step $stride waitfor all
	set nFrames [molinfo $mole get numframes]
	puts [format "Reading %i frames." $nFrames]

	# Move forward, computing
	# current at each step.
	for {set f 0} {$f < $nFrames} {incr f} {
	     molinfo $mole set frame $f

	    # Get the time in nanoseconds for this frame.
	    set t [expr ($nFrames0+$f)*$dt]

	    # Compute the current for each selection.
	    set pot 0.0
	    
	    foreach pInd $probeInd {
		foreach q $charge ind $indList {
		    set d [measure bond [list [list $ind $mole] [list $pInd $probeMol]]]
		    set pot [expr $pot + $q*exp(-$debye*$d)/$d]
		}
	    }
	    # Output in millivolts.
	    set pot [expr $elecConst*$pot/$probeNum*1e10*1000]
	    
	    # Write the total current.
	    puts $out "$t $pot"
	    
	    # Update the display.
	    if {$f % $displayPeriod == 0} {
		puts -nonewline [format "FRAME %i: " $f]
		puts "$t $pot"
	    }
	}
	set nFrames0 [expr $nFrames+$nFrames0]
    }
    close $out
    
    $sel delete
    mol delete $probeMol
    mol delete $mole 
}



