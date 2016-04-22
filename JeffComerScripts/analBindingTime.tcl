# jcomer2@illinois.edu

# jcomer2@illinois.edu
if {$argc < 5} {
    puts "$argv0 name structPrefix outputDir dcdFreq stride dcdFile0 [dcdFile1...]"
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

# Write the binding times.
proc follow {selText bindOn bindOff dt out} {
    set sel [atomselect top $selText]

    # Is the particle initially bound?
    molinfo top set frame 0
    set z [lindex [measure center $sel weight mass] 2]
    if {$z < $bindOn} {
	set isBound 1
	# Skip binding times when we start bound.
	set ready 0
    } else {
	set isBound 0
	set ready 1
    }

    set nFrames [molinfo top get numframes]
    set boundTime 0.0
    for {set f 1} {$f < $nFrames} {incr f} {
	molinfo top set frame $f

	set t [expr {$dt*$f}]
	set z [lindex [measure center $sel weight mass] 2]

	if {$isBound && $z > $bindOff} {
	    # Became unbound.
	    if {$ready} {
		set tau [expr {$t - $boundTime}]
		puts $out "$t $tau"
	    }
	    set ready 1
	    set isBound 0
	} elseif {!$isBound && $z < $bindOn} {
	    # Became bound.
	    set boundTime $t
	    set isBound 1
	}
    }

    $sel delete
    return
}

proc compute {name moiety structPrefix dcdFreq dcdList outDir stride} {
    set bindOn 55
    set bindOff 55
    set displayPeriod 200
    set startFrame 0
    set timestep 1.0
    if {$stride < 1} { set stride 1 }
    set segList "LA65 LA66 LA69 LA70 LA71 LB37 LB39 LB40 LB41 LB42 LB46 LB47 LB48 LB51 LB52 LB53 LB57 LB58 LB59 LB63 LB64 LB65 LB66 LB69 LB70 LB71 LC37 LC39 LC41 LC42 LC46 LC47 LC48 LD51 LD52 LD53 LD55 LD56 LD57 LD58 LD59 LD60 LD61 LD62 LD63 LD64 LD65 LD66 LD67 LD68 LD69 LD70 LD71 LD72 LE37 LE38 LE39 LE40 LE41 LE42 LE43 LE44 LE45 LE46 LE47 LE48 LE49 LE50 LE51 LE52 LE53 LE54 LE55 LE56 LE57 LE58 LE59 LE60 LE61 LE62 LE63 LE64 LE65 LE66 LE67 LE68 LE69 LE70 LE71 LE72 LF37 LF38 LF39 LF40 LF41 LF42 LF43 LF44 LF45 LF46 LF47 LF48 LF49 LF50 LF51 LF52 LF53 LF55 LF56 LF57 LF58 LF64 LG61 LG62 LG67 LG68 LG69 LG72 LH38 LH39 LH43 LH44 LH45 LH46 LH49 LH50 LH51 LH52 LH54 LH55 LH56 LH57 LH60 LH61 LH62 LH63 LH67 LH68 LH69 LH72 LI38 LI39 LI43 LI44 LI45 LI49"
    set selText "name N"
    #set selText "resname PCGL and name N C4 H4A H4B C5 H5A H5B C6 H6A H6B H6C C7 H7A H7B H7C C8 H8A H8B H8C"

    # Input:
    set psf $structPrefix.psf
    set pdb $structPrefix.pdb
    # Output:
    set outPrefix $outDir/time_${name}.dat

    # Get the time change between frames in nanoseconds.
    set dt [expr 1.0e-6*$timestep*$dcdFreq*$stride]

    
    # Load the system.
    mol load psf $psf pdb $pdb
    #set selList {}
    #foreach seg $segList {
#	lappend selList [atomselect top "segname $seg and ($selText)"]
 #   }

    # Loop over the dcd files.
    foreach dcd $dcdList {
	# Load the complete trajectory.
	mol addfile $dcd type dcd step $stride waitfor all
	set n [molinfo top get numframes]
    }
    set nFrames [molinfo top get numframes]
    puts "Read a total of $nFrames frames."
    set totalTime [expr {$nFrames*$dt}]
    puts "Total time: $totalTime"

    set outTotal [open $outDir/total_${name}.dat w]
    puts $outTotal "$totalTime"
    close $outTotal

    # Open the output file.
    set out [open $outPrefix w]

    foreach seg $segList {
	set text "segname $seg and ($selText)"
	follow $text $bindOn $bindOff $dt $out
    }

    close $out

    mol delete top
}

compute $name $moiety $structPrefix $dcdFreq $dcdList $outDir $stride
exit
