# Author: Jeff Comer <jcomer2@illinois.edu>

# Input
#set inTrajList {../8_channel/best10_raw_10mM0.16.traj rest10_raw_10mM1.16.traj}
#set inTrajList {../8_channel/best10_raw_10mM0.16.traj rest10_raw_10mM1.16.traj off10_raw_10mM2.16.traj}
#set inTrajList {../8_channel/best10_middling_10mM0.19.traj rest10_middling_10mM1.19.traj off10_middling_10mM2.19.traj}
#set inTrajList {../8_channel/best10_raw_10mM0.16.traj rest10_raw_10mM1.16.traj}
#set inTrajList {../8_channel/best10_middling_10mM0.19.traj rest10_middling_10mM1.19.traj}

set outDir dcd
set stride 1
set startFrame 0
set nullPos {-1000 0 0}

# Negative frames count from the end.
#set frameList {-1}

proc trimExtension {name} {
    set ind [string last "." $name]
    return [string range $name 0 [expr $ind-1]]
}
proc trimPath {name} {
    set ind [string last "/" $name]
    return [string range $name [expr $ind+1] end]
}

# Construct a pdb line from everything.
proc makePdbLineFull {index segName resId name resName r} {
    set template "ATOM     42  HN  ASN X   4     -41.083  17.391  50.684  0.00  0.00      P1    "

    foreach {x y z} $r {break}
    set record "ATOM  "
    set si [string range [format "     %5i " $index] end-5 end]
    if {[string length $name] < 4} {
	set name [string range " $name    " 0 3]
    } else {
	set name [string range $name 0 3]
    }
    set resName [string range " $resName    " 0 3]
    set temp0 " [string index $segName 0]"
    set resId [string range "    $resId"  end-3 end]
    set temp1 [string range $template  26 29]
    set sx [string range [format "       %8.3f" $x] end-7 end]
    set sy [string range [format "       %8.3f" $y] end-7 end]
    set sz [string range [format "       %8.3f" $z] end-7 end]
    set temp2 [string range $template 54 71]
    set segName [string range "$segName    "  0 3]
    set tempEnd [string range $template 76 end]

    # Construct the pdb line.
    return "${record}${si}${name}${resName}${temp0}${resId}${temp1}${sx}${sy}${sz}${temp2}${segName}${tempEnd}"
}

proc readPositions {fileName} {
    set inp [open $fileName r]
    
    set allFrames {}
    set frame {}
    while {[gets $inp line] >= 0} {
	if {[string match "#*" $line]} {continue}
	if {[string length $line] < 2} {continue}

	if {[string match "END*" $line]} {
	    # End of a frame
	    lappend allFrames $frame
	    set frame {}
	} else {
	    # Append the position data
	    set tok [concat $line]
	    lappend frame [lrange $tok 2 4]
	}
    }

    close $inp
    return $allFrames
}

proc getBiggestFrame {allFrames} {
    set bigLength [llength [lindex $allFrames 0]]
    set bigFrame [lindex $allFrames 0]

    foreach f $allFrames {
	if {[llength $f] > $bigLength} {
	    set bigLength [llength $f] 
	    set bigFrame $f
	}
    }
    return $bigFrame
}

proc writePdb {frame outPdb} {
    # Get the frame.
    set numAtoms [llength $frame]

    set i 1
    set out [open $outPdb w]
    foreach atom $frame {
	foreach {x y z} $atom { break }
	set r [list $x $y $z]
	set atomId [getAtomId $i 0]
	foreach {segName resId name} $atomId { break }

	set line [makePdbLineFull $i $segName $resId $name PAR $r]  
	puts $out $line
	incr i
    }
    close $out

    return
}

proc getAtomId {index type} {
    set nameList {C O N P S}
    
    set nameInd [expr {$type % [llength $nameList]}]
    set name [lindex $nameList $nameInd]
    set resId [expr {$index % 1000}]
    set segNum [expr {$index/1000}]
    
    if {$segNum >= 1000} {
	puts "ERROR: Too many segment names!"
	exit
    }
    set segName "P${segNum}"
    
    return [list $segName $resId $name]
}

proc writeDcd {inTrajList outDir startFrame stride} {
    global nullPos
    # Load the trajectory.
    set allFrames {}
    # Name the resulting file by the last input file.
    set outPrefix [trimExtension $outDir/[trimPath [lindex $inTrajList end]]]
    foreach inTraj $inTrajList {
	set allFrames [concat $allFrames [readPositions $inTraj]]
    }
    puts "Loaded [llength $allFrames] frames."

    # Write a pdb of the biggest frame.
    set bigFrame [getBiggestFrame $allFrames]
    set num [llength $bigFrame]
    writePdb $bigFrame $outPrefix.pdb
    puts "Largest frame has $num atoms."

    # Load this pdb.
    set molId [mol load pdb $outPrefix.pdb]

    set nFrames [llength $allFrames]
    for {set i $startFrame} {$i < $nFrames} {incr i $stride} {
	animate dup $molId
	set frameNum [llength [lindex $allFrames $i]]
	puts "frame $i"

	# Set the particle positions.
	set sel [atomselect top "index < $frameNum"]
	$sel set {x y z} [lindex $allFrames $i]
	$sel delete

	# Set the positions of the leftover particles.
	set sel [atomselect top "index >= $frameNum"]
	$sel set x [lindex $nullPos 0]
	$sel set y [lindex $nullPos 1]
	$sel set z [lindex $nullPos 2]
	$sel delete
    }

    animate delete beg 0 end 0 $molId
    animate write dcd $outPrefix.dcd waitfor all $molId

    puts "Wrote `$outPrefix.dcd'"

    return
}


writeDcd $inTrajList $outDir $startFrame $stride
exit
