# Author: Jeff Comer <jcomer2@illinois.edu>
# Modified by Chen-Yu Li <cli56@illinois.edu>
# Usage: vmd -dispdev text -e $argv0 -args psfPrefix pdbPrefix dcdPrefix xstPrefix dcdFreq timestep selText dimension outPrefix
# Return the current and number of charge carriers in the selection for the frame.
# 2015/8/18


############################################
######Compute current in xyz direction######
############################################

proc computeCurrent {frameCurr sel charge dt} {

    variable minL
    variable maxL
    variable L
    variable dimension

    set frameLast [expr $frameCurr-1]
    
    # Get the current position of the selection.
    $sel frame $frameCurr
    set coor1 [$sel get $dimension]
    
    # Get the last position of the selection.
    $sel frame $frameLast
    set coor0 [$sel get $dimension]
    
    # Get the number of charge carriers.
    set num [$sel num]
    
    # Find the displacements in the z-direction and compute the current.
    set current 0.0
    foreach a0 $coor0 a1 $coor1 {
	
	#Find the absolute value of dZ
	set absdl [expr {abs($a1-$a0)}]

	#Compensate for jumps over the periodic boundary. (Wrap the ion movement)
	set dl [expr $a1-$a0]
	if {$absdl >= 0.5*$L} {
		if {$dl > 0} {
			set dl [expr $dl-$L]
		} elseif {$dl < 0} {
			set dl [expr $dl+$L]
		}
	}

	# Compute the current in nanoamperes.
	set current [expr $current + $charge*$dl/($L*$dt)*1.60217733e-1]
    }
    
    return [list $num $current]
}


proc compute {psfPrefix pdbPrefix dcdPrefix xstPrefix outPrefix} {

    variable minL
    variable maxL
    variable L
    variable dimension
    variable dcdFreq
    variable timestep
    variable selText
    variable R

    set displayPeriod 20
    set stride 1
    set ionText [list "name POT" "name CLA" "name MG"]
    set charge [list 1.0 -1.0 2.0]
    set nameList [list "K" "Cl" "Mg"]

    # Input:
    set psf ${psfPrefix}.psf
    set pdb ${pdbPrefix}.pdb
    set dcd ${dcdPrefix}.dcd
    set xst ${xstPrefix}.xst
    
    # Get the time change between frames in nanoseconds.
    set dt [expr 1.0e-6*$timestep*$dcdFreq*$stride]

    # Open the output files.
    set out {}
    foreach n $nameList {
	lappend out [open "${outPrefix}.curr_${n}.dat" w]
    }
    set outTotal [open "${outPrefix}.curr_total.dat" w]

    # Get the time in nanoseconds.
    set t ""
    set inStream [open ${xstPrefix}.xst r]
    gets $inStream;  gets $inStream;
    while {[gets $inStream line] > 0} {
        lappend t [expr [lindex $line 0] * $timestep * 1.0e-6]
    }

    # Load the system.
    mol load psf $psf pdb $pdb
    #set sel {}
    #foreach st $ionText {
    #    lappend sel [atomselect top $st]
    #}

    set selCenter [atomselect top $selText]
    set all [atomselect top all]


    # Loop over the dcd files.
    # Load the trajectory.
    animate delete all
    mol addfile $dcd type dcd step $stride waitfor all
    set nFrames [molinfo top get numframes]
    puts [format "Reading %i frames." $nFrames]
    
    $selCenter frame 0
    $all frame 0
    
    set r0 [measure center $selCenter weight mass]
    $all moveby [vecinvert $r0]
    
    # Move forward, computing
    # current at each step.
    for {set f 1} {$f < $nFrames} {incr f} {
        $selCenter frame $f
        $all frame $f
    
        set r0 [measure center $selCenter weight mass]
        $all moveby [vecinvert $r0]
    
        # Compute the current for each selection.		
        set currentTotal 0.0
        set currentTotal_in 0.0
        set currentTotal_out 0.0
        foreach ion $ionText q $charge o $out {
    	    #Write the number of carriers and the current for this selection.
            set sel [atomselect top "$ion"]
    	    set data [computeCurrent $f $sel $q $dt]
    	    set current [lindex $data 1]

            set sel_in [atomselect top "$ion and (x*x + y*y < $R*$R)"]
    	    set data_in [computeCurrent $f $sel_in $q $dt]
    	    set current_in [lindex $data_in 1]

            set sel_out [atomselect top "$ion and (x*x + y*y > $R*$R)"]
    	    set data_out [computeCurrent $f $sel_out $q $dt]
    	    set current_out [lindex $data_out 1]

    	    puts $o "[lindex $t $f]\t$current\t$current_in\t$current_out"
    	    
    	    set currentTotal [expr $currentTotal + $current]
    	    set currentTotal_in [expr $currentTotal_in + $current_in]
    	    set currentTotal_out [expr $currentTotal_out + $current_out]

            $sel delete
            $sel_in delete
            $sel_out delete
        }
        
        # Write the total current.
        puts $outTotal "[lindex $t $f]\t$currentTotal\t$currentTotal_in\t$currentTotal_out"
        
        # Update the display.
        if {$f % $displayPeriod == 0} {
    	    puts -nonewline [format "FRAME %i: " $f]
    	    puts "[lindex $t $f]\t$currentTotal"
        }
    }

    foreach o $out {
	close $o
    }
    close $outTotal
    mol delete top
}






#source /home/cli56/scripts/Procs.tcl


if {$argc < 5} {
    puts "vmd -dispdev text -e $argv0 -args psfPrefix pdbPrefix dcdPrefix xstPrefix dcdFreq timestep selText dimension R outPrefix"
    exit
}

set psfPrefix [lindex $argv 0]
set pdbPrefix [lindex $argv 1]
set dcdPrefix [lindex $argv 2]
set xstPrefix [lindex $argv 3]
set dcdFreq [lindex $argv 4]
set timestep [lindex $argv 5]

#map the '_' in selection text to ' '
set selText [lindex $argv 6]
set selText   [string map {_ \ } $selText]

set dimension [lindex $argv 7]
set R [lindex $argv 8]
set outPrefix [lindex $argv 9]


set inStream [open ${xstPrefix}.xst r]
gets $inStream;  gets $inStream;
gets $inStream line
set xs [split $line]

switch -- $dimension {

    "x" {
        set L [lindex $xs 1] 
    }
    "y" {
        set L [lindex $xs 5] 
    }
    "z" {
        set L [lindex $xs 9] 
    }
    default {
        puts "vmd -dispdev text -e $argv0 -args psfPrefix pdbPrefix dcdPrefix xstPrefix dcdFreq timestep selText dimension R outPrefix"
        exit
    }
}

set maxL [expr $L / 2.0]
set minL [expr -1 * $maxL]
puts "$minL\t$maxL"

compute $psfPrefix $pdbPrefix $dcdPrefix $xstPrefix $outPrefix 

exit

