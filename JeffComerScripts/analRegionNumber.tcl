# jcomer2@illinois.edu
if {$argc < 5} {
    puts "$argv0 name moiety structPrefix outputDir dcdFreq stride dcdFile0 \[dcdFile1...\]"
    exit
}
set name [lindex $argv 0]
set moiety [lindex $argv 1]
set structPrefix [lindex $argv 2]
set outDir [lindex $argv 3]
set dcdFreq [lindex $argv 4]
set stride [lindex $argv 5]
set dcdList [lrange $argv 6 end]

source $env(HOME)/scripts/vector.tcl
source $env(HOME)/scripts/gridForce.tcl

proc compute {name moiety structPrefix dcdList dcdFreq outDir stride} {
    set displayPeriod 200
    set timestep 1.0
    set startFrame 0
    set quantList {numK numCl numWater}
    set textList {"name POT" "name CLA" "name OH2"}
    set regionTextList {"z > 42 and z < 75 and x^2+y^2>26^2"}
    set regionNameList {layer}
    if {$stride <= 0} {set stride 1}

    # Open the output files.
    foreach quant $quantList {
	foreach reg $regionNameList {
	    set out($quant,$reg) [open $outDir/${quant}_${name}_${reg}.dat w]
	    puts "$reg $quant"
	}
    }

    # Input:
    set psf $structPrefix.psf
    set pdb $structPrefix.pdb
  
    # Get the time change between frames in nanoseconds.
    set dt [expr 1.0e-6*$timestep*$dcdFreq*$stride]

    # Load the system.
    mol load psf $psf pdb $pdb
    
    # Loop over the dcd files.
    set nFrames0 0
    foreach dcd $dcdList {
	# Load the trajectory.
	animate delete all
	mol addfile $dcd type dcd step $stride waitfor all
	set nFrames [molinfo top get numframes]
	puts [format "Reading %i frames." $nFrames]

	# Move forward computing at every step.
	for {set f $startFrame} {$f < $nFrames} {incr f} {
	    molinfo top set frame $f

	    # Get the time in nanoseconds for this frame.
	    set t [expr {($nFrames0+$f)*$dt}]

	    foreach reg $regionNameList regTex $regionTextList {
		foreach quant $quantList tex $textList {
		    set sel [atomselect top "($tex) and ($regTex)"]
		    set num [$sel num]
		    $sel delete

		    puts $out($quant,$reg) "$t $num"
		}
	    }

	    # Write the time and distance.
	    if {$f % $displayPeriod == 0} {
		puts -nonewline [format "FRAME %i: " $f]
		puts "$t"
	    }
	}
	set nFrames0 [expr $nFrames+$nFrames0]
    }

    foreach quant $quantList {
	foreach reg $regionNameList {
	    close $out($quant,$reg)
	}
    }
    mol delete top
}

compute $name $moiety $structPrefix $dcdList $dcdFreq $outDir $stride
exit
