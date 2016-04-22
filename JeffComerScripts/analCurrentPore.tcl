# jcomer2@illinois.edu
if {$argc < 5} {
    puts "$argv0 name moiety structPrefix outputDir dcdFreq stride dcdFile0 \[dcdFile1...\]"
    exit
}
set name [lindex $argv 0]
set moiety [lindex $argv 1]
set structPrefix [lindex $argv 2]
set outDir [lindex $argv 3]
set dcdFreq [lindex $argv 4]
set startFrame [lindex $argv 5]
set dcdList [lrange $argv 6 end]

source $env(HOME)/scripts/vector.tcl
source $env(HOME)/scripts/gridForce.tcl

 
# Return the current and number of charge carriers in the selection for the frame.
# Does not compensate for jumps over the periodic boundary, although
# they should be ignored for an appropriate choice of cutZ0 and cutZ1.
proc computeCurrent {frameCurr sel charge lz dt cutZ0 cutZ1} {
    set frameLast [expr $frameCurr-1]
    
    # Get the current position of the selection.
    molinfo top set frame $frameCurr
    set z1 [$sel get z]
    
    # Get the last position of the selection.
    molinfo top set frame $frameLast
    set z0 [$sel get z]
    
    # Get the number of charge carriers.
    set num [$sel num]
    
    # Find the displacements in the z-direction and compute the current.
    set currentZ 0.0
    foreach a0 $z0 a1 $z1 {
	# Ignore carriers outside the range for both times.
	if {($a0 < $cutZ0 || $a0 > $cutZ1) && ($a1 < $cutZ0 || $a1 > $cutZ1)} {
	    continue
	}

	# Cut pieces outside of the range.
	if {$a0 < $cutZ0} { set a0 $cutZ0 }
	if {$a1 < $cutZ0} { set a1 $cutZ0 }
	if {$a0 > $cutZ1} { set a0 $cutZ1 }
	if {$a1 > $cutZ1} { set a1 $cutZ1 }
	set dz [expr {$a1-$a0}]

	# Compute the current in nanoamperes.
	set currentZ [expr {$currentZ + $charge*$dz/($lz*$dt)*1.60217733e-1}]
    }
    
    return [list $num $currentZ]
}

proc compute {name moiety structPrefix dcdList dcdFreq outDir startFrame} {
    set displayPeriod 200
    set timestep 1.0
    set displayPeriod 200

    set cutZ0 -20
    set cutZ1 20
    set cutDz [expr $cutZ1-$cutZ0]
    
    set selText [list "name POT" "name CLA"]
    set charge [list 1.0 -1.0]
    set nameList [list "K" "Cl"]
    if {$startFrame < 1} { set startFrame 1 }
    set stride 1
    set pre "thru${cutDz}"

    # Input:
    set psf $structPrefix.psf
    set pdb $structPrefix.pdb
  
    # Get the time change between frames in nanoseconds.
    set dt [expr {1.0e-6*$timestep*$dcdFreq*$stride}]

    # Open the output files.
    set out {}
    foreach n $nameList {
	lappend out [open "${outDir}/${pre}_curr${n}_${name}.dat" w]
    }
    set outTotal [open "${outDir}/${pre}_curr_${name}.dat" w]

    # Load the system.
    mol load psf $psf pdb $pdb
    set sel {}
    foreach st $selText {
	lappend sel [atomselect top $st]
    }

    # Loop over the dcd files.
    set nFrames0 0
    foreach dcd $dcdList {
	# Load the trajectory.
	animate delete all
	mol addfile $dcd type dcd step $stride waitfor all
	set nFrames [molinfo top get numframes]
	puts [format "Reading %i frames." $nFrames]

	# Move forward, computing
	# current at each step.
	for {set f $startFrame} {$f < $nFrames} {incr f} {
	    # Get the time in nanoseconds for this frame.
	    set t [expr ($nFrames0+$f+0.5)*$dt]

	    # Compute the current for each selection.		
	    set currentZTotal 0.0
	    foreach s $sel q $charge o $out {
		#Write the number of carriers and the current for this selection.
		set data [computeCurrent $f $s $q $cutDz $dt $cutZ0 $cutZ1]
		set currentZ [lindex $data 1]
		puts $o "$t $currentZ"
		
		set currentZTotal [expr $currentZTotal + $currentZ]
	    }
	    
	    # Write the total current.
	    puts $outTotal "$t $currentZTotal"
	    
	    # Update the display.
	    if {$f % $displayPeriod == 0} {
		puts -nonewline [format "FRAME %i: " $f]
		puts "$t $currentZTotal"
	    }
	}
	set nFrames0 [expr $nFrames+$nFrames0]
    }


    foreach s $sel o $out {
	$s delete
	close $o
    }
    close $outTotal
    mol delete top
}

compute $name $moiety $structPrefix $dcdList $dcdFreq $outDir $startFrame
exit
