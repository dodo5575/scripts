#This script calculates the number of DNA atoms inside the SiO2 gap.
#Usage: vmd -dispdev text -e $argv0 -args name structPrefix outDir dcdFreq timestep startFrame dcdFile0 \[dcdFile1...\] 
#Written by Chen-Yu Li   cli56@illinois.edu
#Jun. 2 2014



#Compute the number of specific ion in DNA origami at a given frame.
proc computeIon {} {

    set SiO2 [atomselect top "segname U0"]
    set MinMax [measure minmax $SiO2]

    set Zmin [lindex $MinMax 0 2]
    set Zmax [lindex $MinMax 1 2]

    set DNAinPore [atomselect top "nucleic and z>$Zmin and z<$Zmax"]
    
    set num [$DNAinPore num]

    return $num
}

#Main function
proc compute {name structPrefix outDir dcdFreq timestep startFrame dcdList} {

    set displayPeriod 20
    set stride 1

    # Input
    set psf $structPrefix.psf
    set pdb $structPrefix.pdb

    # Get the time change between frames in nanoseconds.
    set dt [expr 1.0e-6*$timestep*$dcdFreq*$stride]

    # Open the output files.
    set out {}
    set outTotal [open "${outDir}/${name}_NumOfDNAatomInPore.dat" w]

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

        set r0 [measure center $nuc weight mass]
        #$all moveby [vecinvert $r0]

        # Move forward, computing
        # current at each step.
        for {set f $startFrame} {$f < $nFrames} {incr f} {

	    animate goto $f

            $nuc frame $f
            $all frame $f

            #set r0 [measure center $nuc weight mass]
            #$all moveby [vecinvert $r0]


            # Get the time in nanoseconds for this frame.
            set t [expr ($nFrames0+$f+0.5)*$dt]

            # Compute the number of ion for each selection.           
            set num [computeIon]
            puts $outTotal "$t $num"


            # Write the total number of ion.

            # Update the display.
            if {$f % $displayPeriod == 0} {
                puts -nonewline [format "FRAME %i: " $f]
                puts "$t $num"
            }
        }
        set nFrames0 [expr $nFrames+$nFrames0]
    }


    close $outTotal
    mol delete top
}

if {$argc < 5} {
    puts "vmd -dispdev text -e $argv0 -args name structPrefix outDir dcdFreq timestep startFrame dcdFile0 \[dcdFile1...\]"
    exit
}

set name [lindex $argv 0]
set structPrefix [lindex $argv 1]
set outDir [lindex $argv 2]
set dcdFreq [lindex $argv 3]
set timestep [lindex $argv 4]
set startFrame [lindex $argv 5]
set dcdList [lrange $argv 6 end]

compute $name $structPrefix $outDir $dcdFreq $timestep $startFrame $dcdList
exit

