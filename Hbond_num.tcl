#This program calculates the distance between two group of atoms, or between two layer of origami.
#Usage: vmd -dispdev text -e $argv0 -args name structPrefix outDir dcdFreq timestep startFrame indexFile limit dcdFile0 \[dcdFile1...\] 
#Attention: use '_' instead of ' ' in selection text
#Chen-Yu Li	cli56@illinois.edu
#2014/12/7

proc computeNum {ind1 ind2 limit f} {

    set num 0
    foreach i1 $ind1 i2 $ind2 {

        set group1 [atomselect top "index $i1"]
        set group2 [atomselect top "index $i2"]

        set v1 [measure center $group1]    
        set v2 [measure center $group2]

        set length [veclength [vecsub $v1 $v2]]
        
        if {$length < $limit} {
            incr num
        }
    }

    return $num
}


proc compute {name structPrefix outDir dcdFreq timestep startFrame indexFile limit dcdList} {

    set displayPeriod 20
    set stride 1

    # Input
    set psf $structPrefix.psf
    set pdb $structPrefix.pdb
    set index [read [open $indexFile r]]

    # parse the index
    set ind1 ""
    set ind2 ""
    for {set i 0} {$i < [llength $index]} {incr i} {
        if {fmod($i,2) == 0} {
            lappend ind1 [lindex $index $i]
        } else {
            lappend ind2 [lindex $index $i]
        }
    }
    puts $ind1
    puts $ind2

    # Get the time change between frames in nanoseconds.
    set dt [expr 1.0e-6*$timestep*$dcdFreq*$stride]

    # Open the output files.
    set out [open "${outDir}/${name}.dat" a+]
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

            set num [computeNum $ind1 $ind2 $limit $f]

            # Write the total number of ion.
            puts $out "$t $num"

            # Update the display.
            if {$f % $displayPeriod == 0} {
                puts -nonewline [format "FRAME %i: " $f]
                puts "$t $num"
            }
        }
        set nFrames0 [expr $nFrames+$nFrames0]
    }

    close $out

    mol delete top
}

if {$argc < 7} {
    puts "vmd -dispdev text -e $argv0 -args name structPrefix outDir dcdFreq timestep startFrame indexFile limit dcdFile0 \[dcdFile1...\]"
    exit
}



set name [lindex $argv 0]
set structPrefix [lindex $argv 1]
set outDir [lindex $argv 2]
set dcdFreq [lindex $argv 3]
set timestep [lindex $argv 4]
set startFrame [lindex $argv 5]
#map the '_' in selection text to ' '
set indexFile [lindex $argv 6]
set limit [lindex $argv 7]
set dcdList [lrange $argv 8 end]

compute $name $structPrefix $outDir $dcdFreq $timestep $startFrame $indexFile $limit $dcdList
exit

