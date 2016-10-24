# Author: Jeff Comer <jcomer2@illinois.edu>
# Modified by Chen-Yu Li <cli56@illinois.edu>
# Usage: vmd -dispdev text -e $argv0 -args psfPrefix pdbPrefix dcd xstPrefix dcdFreq timestep selText dimension outPrefix
# Return the current and number of charge carriers in the selection for the frame.
# 2015/8/18


############################################
######Compute current in xyz direction######
############################################

proc computeCurrent {frameCurr sel charge dt} {

    variable minl
    variable maxl
    variable minL
    variable maxL
    variable l
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
	
        # if both are in the interval
        if {$a0 > $minl && $a0 < $maxl && $a1 > $minl && $a1 < $maxl} {
            set dl [expr $a1-$a0]
        } elseif {$a0 > $minl && $a0 < $maxl} { 
        # if initial coordinate is in the interval
            if {[expr $a1-$a0] > 0} {
                set dl [expr $maxl-$a0]
            } else {
                set dl [expr $minl-$a0]
            }
        } elseif {$a1 > $minl && $a1 < $maxl} {
        # if next coordinate is in the interval
            if {[expr $a1-$a0] > 0} {
                set dl [expr $a1-$minl]
            } else {
                set dl [expr $a1-$maxl]
            }
        } else {
            set dl 0
        }

        # Compute the current in nanoamperes.
        set current [expr $current + $charge*$dl/($l*$dt)*1.60217733e-1]
        
    }
    
    return [list $num $current]
}


proc compute {psfPrefix pdbPrefix dcd xstPrefix outPrefix} {

    variable minl
    variable maxl
    variable minL
    variable maxL
    variable l
    variable L
    variable dimension
    variable dcdFreq
    variable timestep
    variable selText
    variable lipR


    set displayPeriod 20
    set stride 1
    set ionText [list "name POT" "name CLA" "name MG"]
    set charge [list 1.0 -1.0 2.0]
    set nameList [list "K" "Cl" "Mg"]

    # Input:
    set psf ${psfPrefix}.psf
    set pdb ${pdbPrefix}.pdb
    set xst ${xstPrefix}.xst
    
    # Get the time change between frames in nanoseconds.
    set dt [expr 1.0e-6*$timestep*$dcdFreq*$stride]

    # Open the output files.
    set out {}
    foreach n $nameList {
	lappend out [open "${outPrefix}.curr_sec_2_${n}.dat" w]
    }
    set outTotal [open "${outPrefix}.curr_sec_2_total.dat" w]

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

    set selCenter [atomselect top "$selText and same residue as (x*x + y*y < $lipR * $lipR)"]
    set all [atomselect top all]


    # Loop over the dcd files.
    # Load the trajectory.
    animate delete all
    mol addfile $dcd first 0 last -1 waitfor all
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
        foreach selText $ionText q $charge o $out {
    	    #Write the number of carriers and the current for this selection.
            set s [atomselect top "$selText and z > $minL and z < $maxL" frame $f]
    	    set data [computeCurrent $f $s $q $dt]
    	    set current [lindex $data 1]
    	    puts $o "[lindex $t $f]\t$current"
    	    
    	    set currentTotal [expr $currentTotal + $current]
            $s delete
        }
        
        # Write the total current.
        puts $outTotal "[lindex $t $f]\t$currentTotal"
        
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
    puts "vmd -dispdev text -e $argv0 -args psfPrefix pdbPrefix dcd xstPrefix dcdFreq timestep selText dimension l L outPrefix"
    exit
}

set psfPrefix [lindex $argv 0]
set pdbPrefix [lindex $argv 1]
set dcd       [lindex $argv 2]
set xstPrefix [lindex $argv 3]
set dcdFreq   [lindex $argv 4]
set timestep  [lindex $argv 5]

#map the '_' in selection text to ' '
set selText   [lindex $argv 6]
set selText   [string map {_ \ } $selText]

set dimension [lindex $argv 7]
set lipR      [lindex $argv 8]
set l         [lindex $argv 9]
set L         [lindex $argv 10]
set outPrefix [lindex $argv 11]


set inStream [open ${xstPrefix}.xst r]
gets $inStream;  gets $inStream;
gets $inStream line
set xs [split $line]

switch -- $dimension {

    "x" {
        set boxL [lindex $xs 1] 
    }
    "y" {
        set boxL [lindex $xs 5] 
    }
    "z" {
        set boxL [lindex $xs 9] 
    }
    default {
        puts "vmd -dispdev text -e $argv0 -args psfPrefix pdbPrefix dcd xstPrefix dcdFreq timestep selText dimension l L outPrefix"
        exit
    }
}

set maxl [expr $l / 2.0]
set minl [expr -1 * $maxl]
puts "$minl\t$maxl"

set maxL [expr $L / 2.0]
set minL [expr -1 * $maxL]
puts "$minL\t$maxL"

compute $psfPrefix $pdbPrefix $dcd $xstPrefix $outPrefix 

exit

