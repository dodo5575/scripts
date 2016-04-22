#This script calculates the number of ions inside the origami.
#Usage: vmd -dispdev text -e $argv0 -args name structPrefix outDir dcdFreq timestep startFrame dcdFile0 \[dcdFile1...\] 
#Written by Chen-Yu Li   cli56@illinois.edu
#Jun. 21 2013



#Compute the number of specific ion in DNA origami at a given frame.
proc computeIon {selIon} {

    set DNA [atomselect top nucleic]
    set MinMax [measure minmax $DNA]

    set Zmin [lindex $MinMax 0 2]
    set Zmax [lindex $MinMax 1 2]

    if {[expr $Zmax -$Zmin] > [molinfo top get c]} {

	set WT [atomselect top water]

	set zCenter [lindex [measure center $WT] 2]

	set partA [atomselect top "nucleic and z > $zCenter"]
	set partB [atomselect top "nucleic and z < $zCenter"]

	set partAminZ [lindex [measure minmax $partA] 0 2]
	set partBmaxZ [lindex [measure minmax $partB] 1 2]

	set sel [atomselect top "$selIon and (z > $partAminZ or z < $partBmaxZ)"]
    } else {

        set sel [atomselect top "$selIon and z > $Zmin and z < $Zmax"]
    
    }
    
    set num [$sel num]

    return $num
}

#Main function
proc compute {name structPrefix outDir dcdFreq timestep startFrame dcdList} {

    set displayPeriod 20
    set stride 1
    set selText [list "name POT" "name CLA" "name MG"]
    set nameList [list "K" "Cl" "Mg"]

    # Input
    set psf $structPrefix.psf
    set pdb $structPrefix.pdb

    # Get the time change between frames in nanoseconds.
    set dt [expr 1.0e-6*$timestep*$dcdFreq*$stride]

    # Open the output files.
    set out {}
    foreach n $nameList {
        lappend out [open "${outDir}/${name}_NumOf${n}InOrigami.dat" w]
    }
    set outTotal [open "${outDir}/${name}_NumOfIonInOrigami.dat" w]

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
        $all moveby [vecinvert $r0]

        # Move forward, computing
        # current at each step.
        for {set f $startFrame} {$f < $nFrames} {incr f} {

	    animate goto $f

            $nuc frame $f
            $all frame $f

            set r0 [measure center $nuc weight mass]
            $all moveby [vecinvert $r0]


            # Get the time in nanoseconds for this frame.
            set t [expr ($nFrames0+$f+0.5)*$dt]

            # Compute the number of ion for each selection.           
            set TotalIon 0.0
            foreach s $selText o $out {
                #Write the number of carriers and the current for this selection.
                set selIon [computeIon $s]
                puts $o "$t $selIon"

                set TotalIon [expr $TotalIon + $selIon]
            }

            # Write the total number of ion.
            puts $outTotal "$t $TotalIon"

            # Update the display.
            if {$f % $displayPeriod == 0} {
                puts -nonewline [format "FRAME %i: " $f]
                puts "$t $TotalIon"
            }
        }
        set nFrames0 [expr $nFrames+$nFrames0]
    }


    foreach o $out {
	close $o
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

