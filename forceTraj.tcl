# Author: Shu-Han Chao (modified from Jeff's currTraj.tcl)
# Return the force in the selection for the frame.

proc computeforce {k frameCurr sel res direction} {

#    set forceTraj [open "~/schao3/Methylation/analysis/nop/forcetraj/forcez_${frameCurr}.dat" w]

    # get current position of every atom in the selection
    molinfo top set frame $frameCurr
    set position [$sel get $direction]
    set id [$sel get resid]

    set force 0.0
    foreach atom $position atomr $res atomid $id { 
	puts [llength $atom]	
	puts [llength $atomr]	
	set dr [expr $atom - $atomr]
	set force [expr $force + $dr*$k]
#	if {$frameCurr eq 49} {
	    #puts "Frame: $frameCurr"
#	    puts $forceTraj "$atomid $dr"
	    #puts "force: $force"
#	}    
    }

#    close $forceTraj
    return $force
}

proc forcecompute {dcdList} {
    
    # k = 138958 pN/A
    set k 138958 
    set dcdFreq 5000
    set displayPeriod 10
    set stride 1
    set startFrame 1
    set outDir .
    set timestep 2.0
    set seltext "nucleic"
    set direction y
    if {$startFrame < 1} { set startFrame 1 }

    # Input:
    set psf pbcPlate-WI.psf
    set pdb pbcPlate-WI.pdb
    set respdb constrainDNA.pdb
    #set respsf ./sin_${sys}/sin_${sys}_ions.psf
    #set xsc .xsc
    
    # Get the time change between frames in nanoseconds.
    set dt [expr 1.0e-6*$timestep*$dcdFreq*$stride]

    # Open the output files.
    set outTotal [open "~/squarePlate/force${direction}.dat" w]

    # Load the restrain file and read the anchor position.
    mol load pdb $respdb
    set ressel [atomselect top $seltext]
    set res [$ressel get $direction]

    # Load the system.
    mol load psf $psf pdb $pdb
    set sel [atomselect top $seltext]
    
    # Loop over the dcd files.
    set nFrames0 0
    foreach dcdFile $dcdList {
	# Load the trajectory.
	animate delete all
	mol addfile $dcdFile type dcd step $stride waitfor all
	set nFrames [molinfo top get numframes]
	puts [format "Reading %i frames." $nFrames]

	# Move forward, computing
	# Force at each step.
	for {set f $startFrame} {$f < $nFrames} {incr f} {
	    # Get the time in nanoseconds for this frame.
	    set t [expr ($nFrames0+$f+0.5)*$dt]

	    # Compute the Force for each selection.		
	    set forceTotal 0.0
	    set forceTotal [computeforce $k $f $sel $res $direction]
		    
	    # Write the total force in pN.
	    puts $outTotal "$t $forceTotal"
	    
	    # Update the display.
	    if {$f % $displayPeriod == 9} {
		puts -nonewline [format "FRAME %i: " $f]
		puts "$t $forceTotal"
	    }
	}
	set nFrames0 [expr $nFrames+$nFrames0]
    }

    close $outTotal
    mol delete top
}
