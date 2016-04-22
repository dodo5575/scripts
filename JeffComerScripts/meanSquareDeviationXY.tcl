#!/usr/bin/tclsh
# Compute 
# Author: Jeff Comer <jcomer2@illinois.edu>

source $env(HOME)/scripts/vector.tcl

if {$argc < 2} {
    puts "Usage: $argv0 outName timestep inputFile0 inputFile1..."
    exit
}

set outName [lindex $argv 0]
set timestep [lindex $argv 1]
set inFileList [lrange $argv 2 end]

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

proc squareDeviationsXY {posList} {
    set devList {}
    set pos0 [lindex $posList 0]

    foreach pos $posList {
	set d [vecSub $pos $pos0]
	set dev [expr {[lindex $d 0]*[lindex $d 0] + [lindex $d 1]*[lindex $d 1]}]
	lappend devList $dev
    }
    return $devList
}

proc deviationsXY {posList} {
    set devList {}
    set pos0 [lindex $posList 0]

    foreach pos $posList {
	set d [vecSub $pos $pos0]
	set dev [list [lindex $d 0] [lindex $d 1]]
	lappend devList $dev
    }
    return $devList
}

proc squareRoot {valList} {
    set ret {}
    foreach val $valList {
	lappend ret [expr {sqrt($val)}]
    }
    return $ret
}

proc addTimeColumn {valList dt} {
    set ret {}
    set i 0
    foreach val $valList {
	lappend ret [list [expr {$i*$dt}] $val]
	incr i
    }
    return $ret
}

proc differentiate {posList dt} {
    set velList {}
    set n [llength $posList]

    for {set i 1} {$i < $n} {incr i} {
	set v [vecScale [expr {1.0/$dt}] [vecSub [lindex $posList $i] [lindex $posList [expr {$i-1}]]]]
	lappend velList $v
    }
    return $velList
}

proc loadInstances {instancePosVar inFileList} {
    upvar $instancePosVar instancePos

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

	set steps 0
	while {[gets $in line] >= 0} {
	    if {[string length $line] <= 1} { continue }
	    if {[string match "#*" $line]} { continue }
	    
	    # We've reached the end of an instance. Store it.
	    if {[string match "END*" $line]} {
		# Store the position and force lists.
		lappend instancePos $posTraj
		incr count
		
		set posTraj {}
		set forceTraj {}
		set steps 0
		continue
	    }

	    set item [concat $line]
	    if {[llength $item] < 3} { 
		puts "Only [llength $item] of 3 columns found."
		continue 
	    }
	    
	    # Build up the lists.
	    lappend posTraj [lrange $line 0 2]
	    incr steps
	}
	close $in
    }
    return $count
}

# Load the instances (mini-trajectories).
set instPosList {}

set good 0
set inFileListGood {}
foreach inFile $inFileList {
    if {[file exists $inFile]} {
	set good 1
	lappend inFileListGood $inFile
    } else {
	puts "Warning! File `$inFile' does not exist."
    }
}

if {!$good} {
    puts "None of the input files exist!"
    exit
}
puts "Timestep: $timestep"

set nInst [loadInstances instPosList $inFileListGood]
puts "Read $nInst instances."
puts "WEIGHT: $outName $nInst"

############################################
# Compute the square deviations in the xy plane.
set instDevList {}
foreach inst $instPosList {
    lappend instDevList [squareDeviationsXY $inst]
}

# Write the little trajectories.
if {0} {
    set i 0
    foreach inst $instPosList {
	writeList data/$outName.$i.dat [squareDeviationsXY $inst]
	incr i
    }
}

# Get the mean values.
set meanDevList [meanTrajList $instDevList]

# Add the time.
set outList [addTimeColumn $meanDevList $timestep]

# Write the result.
writeList $outName.msd $outList

puts "Complete."
exit
