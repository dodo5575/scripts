# Author: Jeff Comer <jcomer2@illinois.edu>
# Modified by Chen-Yu Li <cli56@illinois.edu>
# Usage: vmd -dispdev text -e $argv0 -args name structPrefix outputDir dcdFreq timestep startFrame dimension minL maxL dcdFile0 \[dcdFile1...\] 
# Return the waterflux in the selection for the frame.


#############################################
######Compute water flux in x direction######
#############################################

proc computeCurrentX {frameCurr sel charge lX dt cutX0 cutX1} {
    set frameLast [expr $frameCurr-1]
    
    # Get the current position of the selection.
    $sel frame $frameCurr
    set X1 [$sel get x]
    
    # Get the last position of the selection.
    $sel frame $frameLast
    set X0 [$sel get x]
    
    # Get the number of charge carriers.
    set num [$sel num]
    
    # Find the displacements in the x-direction and compute the current.
    set currentX 0.0
    foreach a0 $X0 a1 $X1 {
	
	#Find the absolute value of dX
	set absdX [expr {abs($a1-$a0)}]

	#Compensate for jumps over the periodic boundary. (Wrap the ion movement)
	set dX [expr $a1-$a0]
	if {$absdX >= 0.5*$lX} {
		if {$dX > 0} {
			set dX [expr $dX-$lX]
		} elseif {$dX < 0} {
			set dX [expr $dX+$lX]
		}
	}

	# Compute the current in nanoamperes.
	set currentX [expr $currentX + $charge*$dX/($lX*$dt)]
    }
    
    return [list $num $currentX]
}



proc computeX {name structPrefix outDir dcdFreq timestep startFrame cutX0 cutX1 dcdList } {
    set cutDX [expr {$cutX1-$cutX0}]
    set displayPeriod 20
    set stride 1
    set selText [list "resname TIP3 and name OH2"]
    set charge [list 1.0 ]
    set nameList [list "Water"]
    if {$startFrame < 1} { set startFrame 1 }
    

    # Input:
    set psf $structPrefix.psf
    set pdb $structPrefix.pdb
    
    
    # Get the time change between frames in nanoseconds.
    set dt [expr 1.0e-6*$timestep*$dcdFreq*$stride]

    # Open the output files.
    set out {}
    foreach n $nameList {
	lappend out [open "${outDir}/flux${n}_${name}.dat" w]
    }
    

    # Load the system.
    mol load psf $psf pdb $pdb
    set sel {}
    foreach st $selText {
	lappend sel [atomselect top $st]
    }

    set nuc [atomselect top "nucleic"]
    set all [atomselect top all]


    # Loop over the dcd files.
    set nFrames0 0
    foreach dcdFile $dcdList {
	# Load the trajectory.
	animate delete all
	mol addfile $dcdFile type dcd step $stride waitfor all
	set nFrames [molinfo top get numframes]
	puts [format "Reading %i frames." $nFrames]

	$nuc frame 0
        $all frame 0

        set r0 [measure center $nuc weight mass]
        $all moveby [vecinvert $r0]

	# Move forward, computing
	# current at each step.
	for {set f $startFrame} {$f < $nFrames} {incr f} {
            $nuc frame $f
            $all frame $f

            set r0 [measure center $nuc weight mass]
            $all moveby [vecinvert $r0]


	    # Get the time in nanoseconds for this frame.
	    set t [expr ($nFrames0+$f+0.5)*$dt]

	    # Compute the current for each selection.		
	    set currentXTotal 0.0
	    foreach s $sel q $charge o $out {
		#Write the number of carriers and the current for this selection.
		set data [computeCurrentX $f $s $q $cutDX $dt $cutX0 $cutX1]
		set currentX [lindex $data 1]
		puts $o "$t $currentX"
		
		
	    }
	    
	    
	    # Update the display.
	    if {$f % $displayPeriod == 0} {
		puts -nonewline [format "FRAME %i: " $f]
		puts "$t $currentX"
	    }
	}
	set nFrames0 [expr $nFrames+$nFrames0]
    }

    foreach s $sel o $out {
	$s delete
	close $o
    }
    
    mol delete top
}


#############################################
######Compute water flux in y direction######
#############################################

proc computeCurrentY {frameCurr sel charge lY dt cutY0 cutY1} {
    set frameLast [expr $frameCurr-1]
    
    # Get the current position of the selection.
    $sel frame $frameCurr
    set Y1 [$sel get y]
    
    # Get the last position of the selection.
    $sel frame $frameLast
    set Y0 [$sel get y]
    
    # Get the number of charge carriers.
    set num [$sel num]
    
    # Find the displacements in the y-direction and compute the current.
    set currentY 0.0
    foreach a0 $Y0 a1 $Y1 {
	
	#Find the absolute value of dY
	set absdY [expr {abs($a1-$a0)}]

	#Compensate for jumps over the periodic boundary. (Wrap the ion movement)
	set dY [expr $a1-$a0]
	if {$absdY >= 0.5*$lY} {
		if {$dY > 0} {
			set dY [expr $dY-$lY]
		} elseif {$dY < 0} {
			set dY [expr $dY+$lY]
		}
	}

	# Compute the current in nanoamperes.
	set currentY [expr $currentY + $charge*$dY/($lY*$dt)]
    }
    
    return [list $num $currentY]
}



proc computeY {name structPrefix outDir dcdFreq timestep startFrame cutY0 cutY1 dcdList } {
    set cutDY [expr {$cutY1-$cutY0}]
    set displayPeriod 20
    set stride 1
    set selText [list "resname TIP3 and name OH2"]
    set charge [list 1.0 ]
    set nameList [list "Water"]
    if {$startFrame < 1} { set startFrame 1 }
    

    # Input:
    set psf $structPrefix.psf
    set pdb $structPrefix.pdb
    
    
    # Get the time change between frames in nanoseconds.
    set dt [expr 1.0e-6*$timestep*$dcdFreq*$stride]

    # Open the output files.
    set out {}
    foreach n $nameList {
	lappend out [open "${outDir}/flux${n}_${name}.dat" w]
    }
    

    # Load the system.
    mol load psf $psf pdb $pdb
    set sel {}
    foreach st $selText {
	lappend sel [atomselect top $st]
    }

    set nuc [atomselect top "nucleic"]
    set all [atomselect top all]


    # Loop over the dcd files.
    set nFrames0 0
    foreach dcdFile $dcdList {
	# Load the trajectory.
	animate delete all
	mol addfile $dcdFile type dcd step $stride waitfor all
	set nFrames [molinfo top get numframes]
	puts [format "Reading %i frames." $nFrames]

	$nuc frame 0
        $all frame 0

        set r0 [measure center $nuc weight mass]
        $all moveby [vecinvert $r0]

	# Move forward, computing
	# current at each step.
	for {set f $startFrame} {$f < $nFrames} {incr f} {
            $nuc frame $f
            $all frame $f

            set r0 [measure center $nuc weight mass]
            $all moveby [vecinvert $r0]


	    # Get the time in nanoseconds for this frame.
	    set t [expr ($nFrames0+$f+0.5)*$dt]

	    # Compute the current for each selection.		
	    set currentYTotal 0.0
	    foreach s $sel q $charge o $out {
		#Write the number of carriers and the current for this selection.
		set data [computeCurrentY $f $s $q $cutDY $dt $cutY0 $cutY1]
		set currentY [lindex $data 1]
		puts $o "$t $currentY"
		
		
	    }
	    
	    
	    # Update the display.
	    if {$f % $displayPeriod == 0} {
		puts -nonewline [format "FRAME %i: " $f]
		puts "$t $currentY"
	    }
	}
	set nFrames0 [expr $nFrames+$nFrames0]
    }

    foreach s $sel o $out {
	$s delete
	close $o
    }
    
    mol delete top
}


#############################################
######Compute water flux in z direction######
#############################################

proc computeCurrentZ {frameCurr sel charge lZ dt cutZ0 cutZ1} {
    set frameLast [expr $frameCurr-1]
    
    # Get the current position of the selection.
    $sel frame $frameCurr
    set Z1 [$sel get z]
    
    # Get the last position of the selection.
    $sel frame $frameLast
    set Z0 [$sel get z]
    
    # Get the number of charge carriers.
    set num [$sel num]
    
    # Find the displacements in the z-direction and compute the current.
    set currentZ 0.0
    foreach a0 $Z0 a1 $Z1 {
	
	#Find the absolute value of dZ
	set absdZ [expr {abs($a1-$a0)}]

	#Compensate for jumps over the periodic boundary. (Wrap the ion movement)
	set dZ [expr $a1-$a0]
	if {$absdZ >= 0.5*$lZ} {
		if {$dZ > 0} {
			set dZ [expr $dZ-$lZ]
		} elseif {$dZ < 0} {
			set dZ [expr $dZ+$lZ]
		}
	}

	# Compute the current in nanoamperes.
	set currentZ [expr $currentZ + $charge*$dZ/($lZ*$dt)]
    }
    
    return [list $num $currentZ]
}



proc computeZ {name structPrefix outDir dcdFreq timestep startFrame cutZ0 cutZ1 dcdList} {
    set cutDZ [expr {$cutZ1-$cutZ0}]
    set displayPeriod 20
    set stride 1
    set selText [list "resname TIP3 and name OH2"]
    set charge [list 1.0 ]
    set nameList [list "Water"]
    if {$startFrame < 1} { set startFrame 1 }
    

    # Input:
    set psf $structPrefix.psf
    set pdb $structPrefix.pdb
    
    
    # Get the time change between frames in nanoseconds.
    set dt [expr 1.0e-6*$timestep*$dcdFreq*$stride]

    # Open the output files.
    set out {}
    foreach n $nameList {
	lappend out [open "${outDir}/flux${n}_${name}.dat" a+]
    }
    

    # Load the system.
    mol load psf $psf pdb $pdb
    set sel {}
    foreach st $selText {
	lappend sel [atomselect top $st]
    }

    set nuc [atomselect top "nucleic"]
    set all [atomselect top all]


    # Loop over the dcd files.
    set nFrames0 0
    foreach dcdFile $dcdList {
	# Load the trajectory.
	animate delete all
	mol addfile $dcdFile type dcd step $stride waitfor all
	set nFrames [molinfo top get numframes]
	puts [format "Reading %i frames." $nFrames]

	$nuc frame 0
        $all frame 0

        set r0 [measure center $nuc weight mass]
        $all moveby [vecinvert $r0]

	# Move forward, computing
	# current at each step.
	for {set f $startFrame} {$f < $nFrames} {incr f} {
            $nuc frame $f
            $all frame $f

            set r0 [measure center $nuc weight mass]
            $all moveby [vecinvert $r0]


	    # Get the time in nanoseconds for this frame.
	    set t [expr ($nFrames0+$f+0.5)*$dt]

	    # Compute the current for each selection.		
	    set currentZTotal 0.0
	    foreach s $sel q $charge o $out {
		#Write the number of carriers and the current for this selection.
		set data [computeCurrentZ $f $s $q $cutDZ $dt $cutZ0 $cutZ1]
		set currentZ [lindex $data 1]
		puts $o "$t $currentZ"
		
		
	    }
	    
	    
	    # Update the display.
	    if {$f % $displayPeriod == 0} {
		puts -nonewline [format "FRAME %i: " $f]
		puts "$t $currentZ"
	    }
	}
	set nFrames0 [expr $nFrames+$nFrames0]
    }

    foreach s $sel o $out {
	$s delete
	close $o
    }
    
    mol delete top
}



source /home/cli56/scripts/Procs.tcl

if {$argc < 5} {
    puts "vmd -dispdev text -e $argv0 -args name structPrefix outputDir dcdFreq timestep startFrame xscFile dimension dcdFile0 \[dcdFile1...\]"
    exit
}

set name [lindex $argv 0]
set structPrefix [lindex $argv 1]
set outDir [lindex $argv 2]
set dcdFreq [lindex $argv 3]
set timestep [lindex $argv 4]
set startFrame [lindex $argv 5]
set xscFile [lindex $argv 6]
set dimension [lindex $argv 7]
set dcdList [lrange $argv 8 end]

set xs                          [readExtendedSystem $xscFile]

switch -- $dimension {

    "x" {
        set maxL [expr [lindex $xs 1] / 2.0]
        set minL [expr -1 * $maxL]
        puts "$minL\t$maxL"
	computeX $name $structPrefix $outDir $dcdFreq $timestep $startFrame $minL $maxL $dcdList
    }
    "y" {
        set maxL [expr [lindex $xs 5] / 2.0]
        set minL [expr -1 * $maxL]
        puts "$minL\t$maxL"
	computeY $name $structPrefix $outDir $dcdFreq $timestep $startFrame $minL $maxL $dcdList
    }
    "z" {
        set maxL [expr [lindex $xs 9] / 2.0]
        set minL [expr -1 * $maxL]
        puts "$minL\t$maxL"
	computeZ $name $structPrefix $outDir $dcdFreq $timestep $startFrame $minL $maxL $dcdList
    }
    default {
	puts "vmd -dispdev text -e $argv0 -args name structPrefix outDir dcdFreq timestep startFrame xscFile dimension dcdFile0 \[dcdFile1...\]"
    }

}


exit

