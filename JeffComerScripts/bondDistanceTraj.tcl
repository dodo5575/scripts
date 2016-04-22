# Calculate the distance between two CoMs for a trajectory.
# Writes step time(ns) and position (nm) to a file.
# to use: vmd -dispdev text -e stretchDcd.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

proc compute {name moiety structPrefix dcdList dcdFreq outDir startFrame} {
    set displayPeriod 200
    set stride 1
    set timestep 1.0
    set selText "ions"
    # Input:
    set psf $structPrefix.psf
    set pdb $structPrefix.pdb
    set xsc $structPrefix.xsc
    # Output:
    set outPrefix $outDir/dist_${name}.dat

    # Get the time change between frames in nanoseconds.
    set dt [expr 1.0e-6*$timestep*$dcdFreq*$stride]
    
    # Open the output file.
    set out [open $outPrefix w]

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

    # Load the system.
    mol load psf $psf pdb $pdb
    set sel [atomselect top $selText]
    set index0 [lindex [$sel get index] 0]
    set index1 [lindex [$sel get index] 1]
    $sel delete

    set sel0 [atomselect top "index $index0"]
    set sel1 [atomselect top "index $index1"]

    # Loop over the dcd files.
    set nFrames0 1
    foreach dcd $dcdList {
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
	    set dr [vecsub $r1 $r0]
	    foreach {dx dy dz} $dr { break }
	    if {$dx >= 0.5*$lx} {set dx [expr $dx-$lx]}
	    if {$dx < -0.5*$lx} {set dx [expr $dx+$lx]}
	    if {$dy >= 0.5*$ly} {set dy [expr $dy-$ly]}
	    if {$dy < -0.5*$ly} {set dy [expr $dy+$ly]}
	    if {$dz >= 0.5*$lz} {set dz [expr $dz-$lz]}
	    if {$dz < -0.5*$lz} {set dz [expr $dz+$lz]}
	    set d [veclength [list $dx $dy $dz]]

	    # Write the time and position.
	    puts $out "$t $d"
	    
	    if {$f % $displayPeriod == 0} {
		puts -nonewline [format "FRAME %i: " $f]
		puts "$t $d"
	    }
	}
	set nFrames0 [expr $nFrames+$nFrames0]
    }

    $sel0 delete
    $sel1 delete

    mol delete top
    close $out
    return
}
