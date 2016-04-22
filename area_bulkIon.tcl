
#Main function
proc compute {name structPrefix outDir dcdFreq timestep startFrame dcdList} {

    set displayPeriod 20
    set stride 1
    set selText [list "name POT" "name CLA" "name MG" "name P"]
    set nameList [list "K" "Cl" "Mg" "P"]


    # Input
    set psf $structPrefix.psf
    set pdb $structPrefix.pdb

    # Get the time change between frames in nanoseconds.
    set dt [expr 1.0e-6*$timestep*$dcdFreq*$stride]

    # Open the output files.
    set out {}
    foreach n $nameList {
        lappend out [open "${outDir}/${name}_${n}_area_bulkIon.dat" a+]
    }

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

	foreach s $selText o $out {
            for {set f $startFrame} {$f < $nFrames} {incr f} {

                animate goto $f
                $all frame $f

                set x [molinfo top get a]
                set y [molinfo top get b]
                set z [molinfo top get c]
                set upZ [expr $z/2.0 - 5]
                set botZ [expr -1 * $upZ]

                set area [expr $x * $y]
                set bulkV [expr $x * $y * 10]

                set selIon [atomselect top "$s and (z > $upZ or z < $botZ)"]
                set n [$selIon num]
                #1M = 0.000602 molecule/angstrom^3
                set molar [expr $n / ($bulkV * 0.000602)]

                puts $o "$area\t$molar"

            }
	}

        set nFrames0 [expr $nFrames+$nFrames0]
    }

    foreach o $out {
        close $o
    }

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

