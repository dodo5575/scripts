#This program calculate the minmax of specified dimension of DNA origami plate.
#Usage: vmd -dispdev text -e $argv0 -args name structPrefix outDir dcdFreq timestep startFrame dimension dcdFile0 \[dcdFile1...\] 
#Chen-Yu Li	cli56@illinois.edu
#2013/6/26



proc computeHeight {dimension} {

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

    set DNA [atomselect top nucleic]
    set MinMax [measure minmax $DNA]

    set minL [lindex $MinMax 0 $dim]
    set maxL [lindex $MinMax 1 $dim]

    return [expr $maxL - $minL]
}


proc compute {name structPrefix outDir dcdFreq timestep startFrame dimension dcdList} {

    set displayPeriod 20
    set stride 1

    # Input
    set psf $structPrefix.psf
    set pdb $structPrefix.pdb

    # Get the time change between frames in nanoseconds.
    set dt [expr 1.0e-6*$timestep*$dcdFreq*$stride]

    # Open the output files.
    set out [open "${outDir}/${name}_OrigamiMinMax_${dimension}.dat" a+]
#    puts $out "Time OrigamiMinMax"


    # Load the system.
    mol load psf $psf pdb $pdb

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

        # Move forward, computing
        # current at each step.
        for {set f $startFrame} {$f < $nFrames} {incr f} {

            animate goto $f

            $nuc frame $f
            $all frame $f

            # Get the time in nanoseconds for this frame.
            set t [expr ($nFrames0+$f+0.5)*$dt]

	    set height [computeHeight $dimension]

            # Write the total number of ion.
            puts $out "$t $height"

            # Update the display.
            if {$f % $displayPeriod == 0} {
                puts -nonewline [format "FRAME %i: " $f]
                puts "$t $height"
            }
        }
        set nFrames0 [expr $nFrames+$nFrames0]
    }

    close $out

    mol delete top

}

if {$argc < 5} {
    puts "vmd -dispdev text -e $argv0 -args name structPrefix outDir dcdFreq timestep startFrame dimension dcdFile0 \[dcdFile1...\]"
    exit
}

set name [lindex $argv 0]
set structPrefix [lindex $argv 1]
set outDir [lindex $argv 2]
set dcdFreq [lindex $argv 3]
set timestep [lindex $argv 4]
set startFrame [lindex $argv 5]
set dimension [lindex $argv 6]
set dcdList [lrange $argv 7 end]

compute $name $structPrefix $outDir $dcdFreq $timestep $startFrame $dimension $dcdList
exit

