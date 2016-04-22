# Author: Jeff Comer <jcomer2@illinois.edu>

set stride 100
set startFrame 8000
# Input
set inTraj reserve1.0.traj
# Output
set outDir frames

proc trimExtension {name} {
    set ind [string last "." $name]
    return [string range $name 0 [expr $ind-1]]
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

proc readTraj {fileName} {
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
	    lappend frame $tok

	}
    }

    close $inp
    return $allFrames
}

proc makePdbFromFrame {allFrames index outPrefix} {
    set nFrames [llength $allFrames]
    # Negative indices count from the end.
    if {$index < 0} { set index [expr {$nFrames+$index}] }

    if {$index >= $nFrames || $index < -$nFrames} {
	puts "Invalid frame index $index for $nFrames frames."
	exit
    }

    # Get the frame.
    set frame [lindex $allFrames $index]
    puts "TIME [lindex $frame 0 1]"
    set numAtoms [llength $frame]

    set i 1
    set out [open $outPrefix.$index.pdb w]
    foreach atom $frame {
	foreach {type tim x y z} $atom { break }
	set r [list $x $y $z]
	set atomId [getAtomId $i $type]
	foreach {segName resId name} $atomId { break }

	set line [makePdbLineFull $i $segName $resId $name PAR $r]  
	puts $out $line
	incr i
    }
    close $out

    return
}

proc writeFrame {allFrames index outPrefix} {
    set nFrames [llength $allFrames]
    # Negative indices count from the end.
    if {$index < 0} { set index [expr {$nFrames+$index}] }

    if {$index >= $nFrames || $index < -$nFrames} {
	puts "Invalid frame index $index for $nFrames frames."
	exit
    }

    # Get the frame.
    set frame [lindex $allFrames $index]
    puts "TIME [lindex $frame 0 1]"
    set numAtoms [llength $frame]

    set i 1
    set out [open $outPrefix.$index.dat w]
    foreach atom $frame {
	foreach {type tim x y z} $atom { break }
	set r [list $x $y $z]
	puts $out "$x $y $z"
	incr i
    }
    close $out

    return
}

proc writeParticleNumber {allFrames outName} {
    set out [open $outName w]
    foreach f $allFrames {
	set t [lindex $f 0 1]
	set n [llength $f]
	puts $out "$t $n"
    }    
    close $out
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



set prefix $outDir/[trimExtension $inTraj]
set allFrames [readTraj $inTraj]

set nFrames [llength $allFrames]
set frameList {}
for {set i $startFrame} {$i < $nFrames} {incr i $stride} {
    lappend frameList $i
}
set writeFrames [llength $frameList]

foreach frame $frameList {
    puts "Frame $frame"
    writeFrame  $allFrames $frame $prefix 
}

writeParticleNumber $allFrames $prefix.num.dat

exit
