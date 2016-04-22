proc computeDistance {sel dimension} {

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

    set group [atomselect top $sel]

    set l [lindex [measure center $group] $dim]    

    return $l
}


proc compute {name structPrefix outDir dcdFreq timestep startFrame sel dimension dcdList} {

    set displayPeriod 20
    set stride 1

    # Input
    set psf $structPrefix.psf
    set pdb $structPrefix.pdb

    # Get the time change between frames in nanoseconds.
    set dt [expr 1.0e-6*$timestep*$dcdFreq*$stride]

    # Open the output files.
    set out [open "${outDir}/${name}_OrigamiZ.dat" w]
    #puts $out "Time DistanceBetweenLayers"

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

            set height [computeDistance $sel $dimension]

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
    puts "vmd -dispdev text -e $argv0 -args name structPrefix outDir dcdFreq timestep startFrame sel dimension dcdFile0 \[dcdFile1...\]"
    exit
}



set name [lindex $argv 0]
set structPrefix [lindex $argv 1]
set outDir [lindex $argv 2]
set dcdFreq [lindex $argv 3]
set timestep [lindex $argv 4]
set startFrame [lindex $argv 5]
#map the '_' in selection text to ' '
set sel [lindex $argv 6]
set sel [string map {_ \ } $sel]
set dimension [lindex $argv 7]
set dcdList [lrange $argv 8 end]

compute $name $structPrefix $outDir $dcdFreq $timestep $startFrame $sel $dimension $dcdList
exit


