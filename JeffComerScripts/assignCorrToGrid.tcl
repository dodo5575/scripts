#!/usr/bin/tclsh
# Compute 
# Author: Jeff Comer <jcomer2@illinois.edu>

source vector.tcl
source gridForce.tcl

if {$argc != 3} {
    puts "Usage: ./correlationFunctions.tcl gridFile inputFile0 outName"
    exit
}

set gridFile [lindex $argv 0]
set inFileList [list [lindex $argv 1]]
set outName [lindex $argv 2]
set dt 1.0

proc meanListList {listList} {
    set meanList {}
    set n [llength [lindex $listList 0]]
    
    for {set i 0} {$i < $n} {incr i} {
	set sum 0.0
	set count [llength $listList]
	foreach vec $listList {
	    set sum [expr {$sum + [lindex $vec $i]}]
	}
	lappend meanList [expr $sum/$count]
    }
    
    return $meanList
}

proc meanListListVec {listList} {
    set meanList {}
    set n [llength [lindex $listList 0]]
    
    for {set i 0} {$i < $n} {incr i} {
	set sum [vecZero]
	set count [llength $listList]
	foreach vec $listList {
	    set sum [vecadd $sum [lindex $vec $i]]
	}
	lappend meanList [vecScale [expr {1.0/$count}] $sum]
    }
    
    return $meanList
}

proc writeList {fileName l} {
    set out [open $fileName w]
    foreach i $l {
	puts $out $i
    }
    close $out
    return
}

proc writeCorrFunction {fileName l dt} {
    set out [open $fileName w]
    set step 0
    foreach i $l {
	puts $out "$i [expr {$dt*$step}]"
	incr step
    }
    close $out
    return
}

proc loadInstances {instancePosVar instanceForceVar inFileList} {
    upvar $instancePosVar instancePos
    upvar $instanceForceVar instanceForce

    foreach inFile $inFileList {
	# Open the files.
	set in [open $inFile r]
	puts "Read $inFile."

	set posList {}
	set forceList {}

	set count 0
	while {[gets $in line] >= 0} {
	    if {[string length $line] <= 1} { continue }
	    if {[string match "#*" $line]} { continue }
	    
	    # We've reached the end of an instance. Store it.
	    if {[string match "END*" $line]} {
		# Store the position and force lists.
		lappend instancePos $posList
		lappend instanceForce $forceList

		set posList {}
		set forceList {}
		incr count
		continue
	    }

	    set item [concat $line]
	    if {[llength $item] < 7} { 
		puts "Only [llength $item] of 7 columns found."
		continue 
	    }
	    
	    # Build up the lists.
	    lappend posList [lrange $line 1 3]
	    lappend forceList [lrange $line 4 6]
	}
	close $in
    }
    return
}

proc differentiate {posList dt} {
    set velList {}
    set n [llength $velList]

    for {set i 1} {$i < $n} {incr i} {
	set v [vecScale [expr {1.0/$dt}] [vecSub [lindex $posList $i] [lindex $posList [expr $i-1]]]]
	lappend velList $v
    }
    return $velList
}

proc correlationVec {valList} {
    set corr {}
    set v0 [lindex $valList 0]
    foreach v $valList {
	lappend corr [vecDot $v0 $v]
    }
    return $corr
}

# Get the grid geometry.
readDx grid $gridFile
puts "Read the grid with dimensions $grid(nx) $grid(ny) $grid(nz)."

# Load the instances (mini-trajectories).
set instPosList {}
set instForceList {}
loadInstances instPosList instForceList $inFileList
puts "Loaded [llength $instPosList] instances."

# Compute the center for each trajectory.
# Assign it to a place in the grid.
set occupiedGridPoints {}
set i 0
foreach pos $instPosList {
    set cen [vecZero]
    # Compute the mean position.
    foreach r $pos {
	set cen [vecAdd $cen $r]
    }
    set cen [vecScale [expr {1.0/[llength $pos]}] $cen]

    # Use the mean position to set the associated grid point.
    set ind [nearestLatticePosition grid $cen]

    # Make an array entry for each occupied grid point.
    if {[info exists gridInst($ind)]} {
	# The array entry already exists. Append its index.
	lappend gridInst($ind) $i
    } else {
	# This is a new array entry.
	set gridInst($ind) [list $i]
	lappend occupiedGridPoints $ind
    }
    incr i
}
puts $occupiedGridPoints

# Compute the correlation functions for each grid point.
foreach occIndex $occupiedGridPoints {
    set instList $gridInst($occIndex)
    
    # Compute the velocity correlation functions.
    set corrListList {}
    foreach inst $instList {
	# Get the position list from the instance list.
	set posList [lindex $instPosList $inst]
	set velList [differentiate $posList $dt]
	lappend corrListList [correlationVec $velList]
    }

    # Compute the mean velocity correlation function at this grid point.
    set meanCorrVel [meanListList $corrListList]
    writeCorrFunction ${outName}${occIndex}.dat $meanCorrVel $dt
}
puts "Complete."
exit
