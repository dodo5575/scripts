#!/usr/bin/tclsh
# Compute 
# Author: Jeff Comer <jcomer2@illinois.edu>

source vector.tcl
source gridForce.tcl

if {$argc < 1} {
    puts "Usage: ./computeAutocorr.tcl gridFile inputFile0 inputFile1..."
    exit
}

set keepCount 500
set inFileList [lrange $argv 1 end]
set gridFile [lindex $argv 0]
readDx grid $gridFile

proc trimExtension {name} {
    set ind [string last "." $name]
    return [string range $name 0 [expr {$ind-1}]]
}

proc meanTrajList {listList} {
    set meanList {}
    set n [llength [lindex $listList 0]]
    
    for {set i 0} {$i < $n} {incr i} {
	set sum 0.0
	set count [llength $listList]
	foreach s $listList {
	    set sum [expr {$sum + [lindex $s $i]}]
	}
	lappend meanList [expr $sum/$count]
    }
    
    return $meanList
}

proc meanTrajListVec {listList} {
    set meanList {}
    set n [llength [lindex $listList 0]]
    
    for {set i 0} {$i < $n} {incr i} {
	set sum [vecZero]
	set count [llength $listList]
	foreach vec $listList {
	    set sum [vecAdd $sum [lindex $vec $i]]
	}
	lappend meanList [vecScale [expr {1.0/$count}] $sum]
    }
    
    return $meanList
}

proc addTrajListVec {listList vectorList scale} {
    set retTrajList {}
    foreach lis $listList {
	set retList {}
	foreach vector $lis vector0 $vectorList {
	    set v [vecAdd $vector [vecScale $scale $vector0]]
	    lappend retList [vecAdd $vector [vecScale $scale $vector0]]
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

    global keepCount

    foreach inFile $inFileList {
	# Open the files.
	set in [open $inFile r]
	#puts "Read $inFile."

	set posTraj {}
	set forceTraj {}

	set count 0
	set steps 0
	while {[gets $in line] >= 0} {
	    if {[string length $line] <= 1} { continue }
	    if {[string match "#*" $line]} { continue }
	    
	    # We've reached the end of an instance. Store it.
	    if {[string match "END*" $line]} {
		if {$steps == $keepCount} {
		    # Store the position and force lists.
		    lappend instancePos $posTraj
		    lappend instanceForce $forceTraj
		    incr count
		}

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
    return
}

# Load the instances (mini-trajectories).
foreach inFile $inFileList {
    set instPosList {}
    set instForceList {}

    set dt [readTimestep $inFile]
    loadInstances instPosList instForceList [list $inFile]
    puts "$inFile [llength $instPosList] instances"    

    # Compute the velocity.
    set instVelList {}
    foreach posTraj $instPosList {
	lappend instVelList [differentiate [wrapToSelf grid $posTraj] $dt]
    }

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
    set outName "[trimExtension $inFile].vcorr"
    writeCorrFunction $outName $meanVelCorrFunc $dt

    ##########################################
    # Compute the force correlation functions.
    # Get the mean force over the trajectories <F(r,t)>.
    set meanForceTraj [meanTrajListVec $instForceList]

    # Compute the deviation of the force F(r,t) - <F(r,t)>.
    set deltaForceTrajList [addTrajListVec $instForceList $meanForceTraj -1.0]
    	
    # Compute the force correlation functions.
    set forceCorrFuncList {}
    foreach deltaForceTraj $deltaForceTrajList {
	set correlationFunc [correlationVec $deltaForceTraj]
	lappend forceCorrFuncList $correlationFunc
    }

    # Compute the mean force correlation function at this grid point.
    set meanForceCorrFunc [meanTrajList $forceCorrFuncList]

    # Write the result.
    set outName "[trimExtension $inFile].fcorr"
    writeCorrFunction $outName $meanForceCorrFunc $dt
}

puts "Complete."
exit
