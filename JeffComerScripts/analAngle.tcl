# jcomer2@illinois.edu

# jcomer2@illinois.edu
if {$argc < 5} {
    puts "$argv0 name structPrefix outputDir stride dcdFile0 [dcdFile1...]"
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

proc compute {name moiety structPrefix dcdFreq dcdList outDir stride} {
    set displayPeriod 200
    set radMin 26
    set startFrame 0
    set timestep 1.0
    if {$stride < 1} { set stride 1 }
    set segList "LA65 LA66 LA69 LA70 LA71 LB37 LB39 LB40 LB41 LB42 LB46 LB47 LB48 LB51 LB52 LB53 LB57 LB58 LB59 LB63 LB64 LB65 LB66 LB69 LB70 LB71 LC37 LC39 LC41 LC42 LC46 LC47 LC48 LD51 LD52 LD53 LD55 LD56 LD57 LD58 LD59 LD60 LD61 LD62 LD63 LD64 LD65 LD66 LD67 LD68 LD69 LD70 LD71 LD72 LE37 LE38 LE39 LE40 LE41 LE42 LE43 LE44 LE45 LE46 LE47 LE48 LE49 LE50 LE51 LE52 LE53 LE54 LE55 LE56 LE57 LE58 LE59 LE60 LE61 LE62 LE63 LE64 LE65 LE66 LE67 LE68 LE69 LE70 LE71 LE72 LF37 LF38 LF39 LF40 LF41 LF42 LF43 LF44 LF45 LF46 LF47 LF48 LF49 LF50 LF51 LF52 LF53 LF55 LF56 LF57 LF58 LF64 LG61 LG62 LG67 LG68 LG69 LG72 LH38 LH39 LH43 LH44 LH45 LH46 LH49 LH50 LH51 LH52 LH54 LH55 LH56 LH57 LH60 LH61 LH62 LH63 LH67 LH68 LH69 LH72 LI38 LI39 LI43 LI44 LI45 LI49"
    set selText0 "name N"
    set selText1 "name P"
    
    # Input:
    set psf $structPrefix.psf
    set pdb $structPrefix.pdb
    # Output:
    set outPrefix $outDir/angle_${name}.dat

    # Get the time change between frames in nanoseconds.
    set dt [expr {1.0e-6*$timestep*$dcdFreq*$stride}]
    set radMinSq [expr {$radMin*$radMin}]
    set degreeConst [expr {180.0/(4.0*atan(1.0))}]

    # Open the output file.
    set out [open $outPrefix w]

    # Load the system.
    mol load psf $psf pdb $pdb

    # Make the selections.
    set selList0 {}
    set selList1 {}
    foreach seg $segList {
	lappend selList0 [atomselect top "segname $seg and ($selText0)"]
	lappend selList1 [atomselect top "segname $seg and ($selText1)"]
    }
    
    # Loop over the dcd files.
    set nFrames0 0
    set angle -1.0
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

	    foreach s0 $selList0 s1 $selList1 {
		set r0 [lindex [$s0 get {x y z}] 0]
		set r1 [lindex [$s1 get {x y z}] 0]
		
		foreach {x0 y0 z0} $r0 { break }
		# Check that the lipid is not in the pore opening.
		if {$x0*$x0 + $y0*$y0 > $radMinSq} {
		    # Compute the angle.
		    set d [vecsub $r1 $r0]
		    set dMag [veclength $d]
		    set dUnit [vecscale [expr {1.0/$dMag}] $d]
		    set angle [expr {$degreeConst*acos([lindex $dUnit 2])}]
		    puts $out "$z0 $angle"
		}
	    }
	    
	    if {$f % $displayPeriod == 0} {
		puts [format "FRAME %i: " $f $angle]
	    }
	}
	set nFrames0 [expr {$nFrames+$nFrames0}]
    }
    close $out

    foreach s0 $selList0 s1 $selList1 {
	$s0 delete
	$s1 delete
    }
    mol delete top
}

compute $name $moiety $structPrefix $dcdFreq $dcdList $outDir $stride
exit
