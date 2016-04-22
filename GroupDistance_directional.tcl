#This program calculates the distance between two group of atoms, or between two layer of origami.
#Usage: vmd -dispdev text -e $argv0 -args name structPrefix outDir xstPrefix timestep sel1 sel2 dcdFile0 \[dcdFile1...\] 
#Attention: use '_' instead of ' ' in selection text
#Chen-Yu Li	cli56@illinois.edu
#2015/4/1

proc computeDistance {sel1 sel2} {

    variable dim

    set group1 [atomselect top $sel1]
    set group2 [atomselect top $sel2]

    set l1 [lindex [measure center $group1] $dim]    
    set l2 [lindex [measure center $group2] $dim]

    set length [expr abs($l1 - $l2)]

    return $length
}


proc compute {psfPrefix pdbPrefix dcdPrefix xstPrefix sel1 sel2 timestep outPut} {

    set displayPeriod 20
    set stride 1

    # Input
    set psf $psfPrefix.psf
    set pdb $pdbPrefix.pdb
    set dcd $dcdPrefix.dcd

    # Get the time change between frames in nanoseconds.
    set t ""
    set inStream [open ${xstPrefix}.xst r]
    gets $inStream;  gets $inStream;
    while {[gets $inStream line] > 0} {
        lappend t [expr [lindex $line 0] * $timestep * 1.0e-6]
    }

    # Open the output files.
    set out [open "${outPut}" a+]
    #puts $out "Time DistanceBetweenLayers"

    # Load the system.
    mol load psf $psf pdb $pdb

    set nuc [atomselect top "nucleic"]
    set all [atomselect top all]

    # Loop over the dcd files.
    set nFrames0 0
    foreach dcdFile $dcd {
        # Load the trajectory.
        animate delete all
        mol addfile $dcdFile type dcd step $stride waitfor all
        set nFrames [molinfo top get numframes]
        puts [format "Reading %i frames." $nFrames]

        $nuc frame 0
        $all frame 0

        # Move forward, computing
        # current at each step.
        for {set f 0} {$f < $nFrames} {incr f} {

            animate goto $f

            $nuc frame $f
            $all frame $f

            set height [computeDistance $sel1 $sel2]

            # Write the total number of ion.
            puts $out "[lindex $t $f]\t$height"

            # Update the display.
            if {$f % $displayPeriod == 0} {
                puts -nonewline [format "FRAME %i: " $f]
                puts "[lindex $t $f]\t$height"
            }
        }
        set nFrames0 [expr $nFrames+$nFrames0]
    }

    close $out

    mol delete top
}

if {$argc < 8} {
    puts "vmd -dispdev text -e $argv0 -args psfPrefix pdbPrefix dcdPrefix xstPrefix sel1 sel2 timestep dimension outPut"
    exit
}


set psfPrefix [lindex $argv 0]
set pdbPrefix [lindex $argv 1]
set dcdPrefix [lindex $argv 2]
set xstPrefix [lindex $argv 3]
#map the '_' in selection text to ' '
set sel1 [lindex $argv 4]
set sel1 [string map {_ \ } $sel1]
set sel2 [lindex $argv 5]
set sel2 [string map {_ \ } $sel2]

set timestep [lindex $argv 6]

set dimension [lindex $argv 7]
# convert dimension
variable dim
switch -- $dimension {
    "x" {
    set dim 0
    }
    "y" {
    set dim 1
    }
    "z" {
    set dim 2
    }
    default {
    }
}

set outPut [lindex $argv 8]

compute $psfPrefix $pdbPrefix $dcdPrefix $xstPrefix $sel1 $sel2 $timestep $outPut
exit

