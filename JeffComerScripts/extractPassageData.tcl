#!/usr/bin/tclsh
# Author: Jeff Comer <jcomer2@illinois.edu>

if {$argc < 1} {
    puts "Usage: $argv0 inputFile0 \[inputFile1\]..."
    exit
}

set inFileList [lrange $argv 0 end]

proc trimExtension {name} {
    set ind [string last "." $name]
    return [string range $name 0 [expr {$ind-1}]]
}

proc trimPath {name} {
    set ind [string last "/" $name]
    return [string range $name [expr $ind+1] end]
}

foreach inFile $inFileList {
    # Open the files.
    set in [open $inFile r]
    set outName [trimPath $inFile]

    set outBound [open bound_${outName} w]
    set outFree [open free_${outName} w]
    set outEvent [open event_${outName} w]

    while {[gets $in line] >= 0} {
	if {[string length $line] <= 1} {continue}
	if {[string match "#*" $line]} {continue}

	set tok [concat $line]
	
	set id [lindex $tok 0]
	set enterTime [lindex [split [lindex $tok 1] ","] 0]
	set firstBoundTime [lindex [split [lindex $tok 2] ","] 0]; # This is exitTime if no bound events occurred
	set exitTime [lindex [split [lindex $tok end-1] ","] 0]
 	set passageTime [lindex $tok end]
	set eventTimeList [lrange $tok 2 end-2]
	
	set eventNum [llength $eventTimeList]
	set boundTimeList $eventTimeList
	# Ignore the first event, which marks a binding.
	set freeTimeList [lrange $eventTimeList 1 end]

	# Was the particle bound on entry?
	if {$firstBoundTime == $enterTime} {
	    set enterBound 1
	    # Ignore the first binding event in calculating statistics.
	    set boundTimeList [lrange $boundTimeList 2 end]
	} else {
	    set enterBound 0
	}

	# Was the particle bound on exit?
	if {$eventNum % 2 == 1} {
	    set exitBound 1
	    # Ignore the last bound event that does not finish.
	    set boundTimeList [lrange $boundTimeList 0 end-1]
	} else {
	    # Ignore the last free event that does not finish.
	    set exitBound 0
	    set freeTimeList [lrange $freeTimeList 0 end-1]
	}

	# How many binding and unbinding events were recorded?
	if {$eventNum == 0} {
	    set boundNum 0
	    set freeNum 0
	} elseif {$eventNum == 1 && $enterBound} {
	    set boundNum 0
	    set freeNum 0
	} elseif {$eventNum == 1 && !$enterBound} {
	    set boundNum 1
	    set freeNum 0
	} elseif {$enterBound} {
	    set boundNum [expr {($eventNum+1)/2 - $enterBound}]
	    set freeNum [expr {$eventNum/2}]
	}
     
	if {[llength $boundTimeList] % 2 == 1 || [llength $freeTimeList] % 2 == 1} {
	    puts "Bug in forming bound and free event lists!"
	    exit
	}

	# Write the bound events.
	# Format: "unbindTime boundDuration"
	set boundTSum 0.0
	set boundXSum 0.0
	foreach {e0 e1} $boundTimeList {
	    foreach {t0 x0} [split $e0 ","] { break }
	    foreach {t1 x1} [split $e1 ","] { break }

	    set dt [expr {$t1-$t0}]
    	    set dx [expr {$x1-$x0}]
	    set boundTSum [expr {$boundTSum + $dt}]
	    set boundXSum [expr {$boundXSum + $dx}]
	    puts $outBound "$t1 $dt $dx"
	}
	set boundMeanTime [expr {$boundTSum/([llength $boundTimeList]/2)}]
	set boundMeanDisp [expr {$boundXSum/([llength $boundTimeList]/2)}]

	# Write the free events.
	# Format: "rebindTime freeDuration"
	set freeTSum 0.0
	set freeXSum 0.0
	foreach {e0 e1} $freeTimeList {
	    foreach {t0 x0} [split $e0 ","] { break }
	    foreach {t1 x1} [split $e1 ","] { break }

	    set dt [expr {$t1-$t0}]
    	    set dx [expr {$x1-$x0}]
	    set freeTSum [expr {$freeTSum + $dt}]
	    set freeXSum [expr {$freeXSum + $dx}]
	    puts $outFree "$t1 $dt $dx"
	}
	set freeMeanTime [expr {$freeTSum/([llength $freeTimeList]/2)}]
	set freeMeanDisp [expr {$freeXSum/([llength $freeTimeList]/2)}]

	set meanTime [expr {0.5*($enterTime + $exitTime)}]

	# Write the passage statistics.
	puts $outEvent "$meanTime $enterBound $exitBound $passageTime $boundMeanTime $boundMeanDisp $boundNum $freeMeanTime $freeMeanDisp $freeNum"
	#puts$outEvent " meanTime1 enterBound2 exitBound3 passageTime4 boundMeanTime5 boundMeanDisp6 boundNum7 freeMeanTime8 freeMeanDisp9 freeNum10"
    }

    close $in
    close $outBound
    close $outFree
    close $outEvent
}
