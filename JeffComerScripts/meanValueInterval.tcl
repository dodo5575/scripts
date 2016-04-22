#!/usr/bin/tclsh
# Author: Jeff Comer <jcomer2@illinois.edu>

if {$argc < 2} {
    puts "Usage: tclsh meanValueCut cutTime0 cutTime1 inputFile0 \[inputFile1\]..."
    exit
}

# Input:
set cutTime0 [lindex $argv 0]
set cutTime1 [lindex $argv 1]
set inFileList [lrange $argv 2 end]

if {![string is double $cutTime0]} {
    puts "ERROR: cutTime0 must be a real number."
    exit
}


if {![string is double $cutTime1]} {
    puts "ERROR: cutTime1 must be a real number."
    exit
}


set dimList {val}
foreach inFile $inFileList {
    # Open the files.
    set in [open $inFile r]
    
    foreach dim $dimList {
	set sum($dim) 0
	set sumSq($dim) 0
    }
    set n 0

    while {[gets $in line] >= 0} {
	if {[string length $line] <= 1} {break}
	
	if {[string match "#*" $line]} {continue}
	set tok [concat $line]
	set x(tim) [expr double([lindex $tok 0])]
	set x(val) [expr double([lindex $tok 1])]

	if {$x(tim) >= $cutTime0 && $x(tim) < $cutTime1} {
	    foreach dim $dimList {
		set sum($dim) [expr $sum($dim) + $x($dim)]
		set sumSq($dim) [expr $sumSq($dim) + $x($dim)*$x($dim)]
	    }
	    incr n
	}
    }
    
    #puts "\nRead $inFile"
    #puts "final time: $x(tim)"
    #puts "total time: [expr {$x(tim)-$cutTime}]"

    if {$n == 0} {
	puts "No data."
    } else {
	foreach dim $dimList {
	    set mean($dim) [expr $sum($dim)/$n]
	    set err($dim) [expr sqrt(($sumSq($dim) - $sum($dim)*$sum($dim)/$n)/($n-1)/$n)]
	    #puts "mean $dim: $mean($dim) +/- $err($dim)"
	    puts "$mean($dim) $err($dim)"
	}
    }
    
    close $in
}
