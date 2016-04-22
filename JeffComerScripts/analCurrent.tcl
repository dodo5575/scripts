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
    if {$stride <= 0} {set stride 1}

    # Input:
    set psf $structPrefix.psf
    set xsc electric1.restart.xsc 
  
    # Get the time change between frames in nanoseconds.
    set dt [expr {1.0e-6*$timestep*$dcdFreq*$stride}]

    # Load the system.
    mol load psf $psf
    set sel [atomselect top "ions"]
    set chargeList [$sel get charge]
    
    newGridXsc grid $xsc 3.0
    
    # Calculate the coupling.
    set lz [expr {$grid(nz)*[lindex $grid(delta) 2 2]}]
    set currFactorList {}
    foreach charge $chargeList {
	lappend currFactorList [expr {$charge/($lz*$dt)*1.60217733e-1}]
    }

    set out [open $outDir/curr_${name}.dat w]

    # Loop over the dcd files.
    set nFrames0 0
    foreach dcd $dcdList {
	# Load the trajectory.
	animate delete all
	mol addfile $dcd type dcd step $stride waitfor all
	set nFrames [molinfo top get numframes]
	puts [format "Reading %i frames." $nFrames]

	# Move forward computing at every step.
	for {set f 1} {$f < $nFrames} {incr f} {
	    molinfo top set frame [expr {$f-1}]
	    set pos0List [$sel get {x y z}]
	    molinfo top set frame $f
    	    set posList [$sel get {x y z}]

	    set current 0.0
	    foreach pos $posList pos0 $pos0List q $currFactorList {
		set dz [lindex [wrapDiff grid [vecsub $pos $pos0]] 2]
		set current [expr {$current + $q*$dz}]
	    }

	    # Get the time in nanoseconds for this frame.
	    set t [expr {($nFrames0+$f)*$dt}]
	    puts $out "$t $current"
	    
	    # Write the time and distance.
	    if {$f % $displayPeriod == 0} {
		puts -nonewline [format "FRAME %i: " $f]
		puts "$t $current"
	    }
	}
	set nFrames0 [expr $nFrames+$nFrames0]
    }

    $sel delete
    mol delete top
    close $out
}

compute $name $moiety $structPrefix $dcdList $dcdFreq $outDir $stride
exit
