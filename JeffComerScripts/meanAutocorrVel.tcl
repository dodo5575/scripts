#!/usr/bin/tclsh
# Compute 
# Author: Jeff Comer <jcomer2@illinois.edu>

source $env(HOME)/scripts/vector.tcl

if {$argc < 2} {
    puts "Usage: ./meanAutocorrVel.tcl outName inputFile0 inputFile1..."
    exit
}

set outName [lindex $argv 0]
set inFileList [lrange $argv 1 end]

proc trimExtension {name} {
    set ind [string last "." $name]
    return [string range $name 0 [expr {$ind-1}]]
}

proc meanTrajList {listList} {
    set meanList {}
    set n [llength [lindex $listList 0]]
    set goodListList {}

    foreach lis $listList {
	if {[llength $lis] == $n} {
	    lappend goodListList $lis
	}
    }

    for {set i 0} {$i < $n} {incr i} {
	set sum 0.0
	set count [llength $goodListList]
	foreach s $goodListList {
	    set sum [expr {$sum + [lindex $s $i]}]
	}
	if {$count == 0} {
	    lappend meanList 0.0
	} else {
	    lappend meanList [expr {$sum/$count}]
	}
	
    }
    
    return $meanList
}

proc meanTrajListVec {listList} {
    set meanList {}
    set n [llength [lindex $listList 0]]

    set goodListList {}
    foreach lis $listList {
	if {[llength $lis] == $n} {
	    lappend goodListList $lis
	}
    }

    for {set i 0} {$i < $n} {incr i} {
	set sum [vecZero]
	set count [llength $goodListList]
	foreach vec $goodListList {
	    set sum [vecAdd $sum [lindex $vec $i]]
	}
	if {$count == 0} {
	    lappend meanList [vecZero]
	} else {
	    lappend meanList [vecScale [expr {1.0/$count}] $sum]
	}
    }
    
    return $meanList
}

proc addTrajListVec {listList vectorList scale} {
    set retTrajList {}
    foreach lis $listList {
	set retList {}
	foreach vector $lis vector0 $vectorList {
	    if {[llength $vector] == 3 && [llength $vector0] == 3} {
		set v [vecAdd $vector [vecScale $scale $vector0]]
		lappend retList [vecAdd $vector [vecScale $scale $vector0]]
	    }
	}
	lappend retTrajList $retList
    }
    return $retTrajList
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
	puts $out "[expr {$dt*$step}] $i"
	incr step
    }
    close $out
    return
}

proc differentiate {posList dt} {
    set velList {}
    set n [llength $posList]

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

proc wrapToSelf {gridVar posList} {
    upvar $gridVar grid

    if {[llength $posList] == 0} { return $posList }
    set r0 [lindex $posList 0]

    set ret {}
    foreach pos $posList {
	set dr [vecSub $pos $r0]
	set r1 [vecAdd $r0 [wrapDisplacement grid $dr]]
	lappend ret $r1
    }
    return $ret
}

proc readTimestep {inFile} {
    # Open the files.
    set in [open $inFile r]
    set count 0
    set timList {}
    while {[gets $in line] >= 0 && $count < 2} {
	if {[string match "#*" $line]} { continue }

	set item [concat $line]
	if {[llength $item] < 7} { 
	    puts "Only [llength $item] of 7 columns found."
	    continue 
	}
	lappend timList [lindex $item 0]
	incr count
    }
    close $in

    return [expr {[lindex $timList 1]-[lindex $timList 0]}]
}

proc loadInstances {instancePosVar instanceForceVar inFileList} {
    upvar $instancePosVar instancePos
    upvar $instanceForceVar instanceForce

    set count 0
    foreach inFile $inFileList {
	if {![file exists $inFile]} { 
	    #puts "Warning! $inFile does not exist."
	    continue
	}
	
	# Open the files.
	set in [open $inFile r]
	puts "Read $inFile."

	set posTraj {}
	set forceTraj {}

	set steps 0
	while {[gets $in line] >= 0} {
	    if {[string length $line] <= 1} { continue }
	    if {[string match "#*" $line]} { continue }
	    
	    # We've reached the end of an instance. Store it.
	    if {[string match "END*" $line]} {
		# Store the position and force lists.
		lappend instancePos $posTraj
		#lappend instanceForce $forceTraj
		incr count
		

		set posTraj {}
		set forceTraj {}
		set steps 0
		continue
	    }

	    set item [concat $line]
	    if {[llength $item] < 7} { 
		puts "Only [llength $item] of 7 columns found."
		continue 
	    }
	    
	    # Build up the lists.
	    lappend posTraj [lrange $line 1 3]
	    lappend forceTraj [lrange $line 4 6]
	    incr steps
	}
	close $in
    }
    return $count
}

# Load the instances (mini-trajectories).
set instVelList {}
set instForceList {}

set good 0
set inFileListGood {}
foreach inFile $inFileList {
    if {[file exists $inFile]} {
	set good 1
	set dt [readTimestep $inFile]
	lappend inFileListGood $inFile
    } else {
	puts "Warning! File `$inFile' does not exist."
    }
}

if {!$good} {
    puts "None of the input files exist!"
    exit
}
puts "Timestep: $dt"


set nInst [loadInstances instVelList instForceList $inFileListGood]
puts "Read $nInst instances."
puts "WEIGHT: $outName $nInst"

############################################
# Compute the velocity correlation functions.
# Get the mean velocity over the trajectories <v(r,t)>.
set meanVelTraj [meanTrajListVec $instVelList]

# Compute the deviation of the velocity v(r,t) - <v(r,t)>.
set deltaVelTrajList [addTrajListVec $instVelList $meanVelTraj -1.0]

# Compute the vel correlation functions.
set velCorrFuncList {}
foreach deltaVelTraj $deltaVelTrajList {
    set correlationFunc [correlationVec $deltaVelTraj]
    lappend velCorrFuncList $correlationFunc
}

# Compute the mean velocity correlation function at this grid point.
set meanVelCorrFunc [meanTrajList $velCorrFuncList]

# Write the result.
writeCorrFunction $outName.vcorr $meanVelCorrFunc $dt

puts "Complete."
exit
