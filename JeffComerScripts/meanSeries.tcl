#!/usr/bin/tclsh
# Author: Jeff Comer <jcomer2@illinois.edu>

set fileList $argv
set outPrefix mean

proc trimExtension {name} {
    set ind [string last "." $name]
    return [string range $name 0 [expr $ind-1]]
}
proc trimPath {name} {
    set ind [string last "/" $name]
    return [string range $name [expr $ind+1] end]
}
proc trimPrefix {name} {
    set ind [string last "_" $name]
    return [string range $name [expr $ind+1] end]
}
# Just read a space delimited data file.
proc readData {fileName} {
    set in [open $fileName r]
    
    set r {}
    while {[gets $in line] >= 0} {
	if {[string match "#*" $line]} {continue}
	if {[string length $line] < 2} {continue}

	set tok [concat $line]
	lappend r [lrange $tok 0 1]
    }

    close $in
    return $r
}


# Read the files.
set trajList {}
foreach file $fileList {
    lappend trajList [readData $file]
}
set trajNum [llength $trajList]
puts "Found $trajNum files for this set."

# Find the longest trajectory.
set long [llength [lindex $trajList 0]]
foreach traj $trajList {
    if {[llength $traj] > $long} { set long [llength $traj] }
}
puts "Longest trajectory: $long frames"

# Open the output file.
set out [open ${outPrefix}_[lindex $fileList end] w]

# Compute the average.
# Run over each frame.
for {set i 0} {$i < $long} {incr i} {
    set sum 0.0
    set sumSq 0.0
    set n 0
    # Run over each trajectory.
    for {set j 0} {$j < $trajNum} {incr j} {
	set val [lindex $trajList $j $i 1]
	if {[llength $val] > 0} {
	    set tim [lindex $trajList $j $i 0]
	    set sum [expr {$sum + $val}]
	    set sumSq [expr {$sumSq + $val*$val}]
	    incr n
	}
    }

    # Compute the mean and error.
    set mean [expr {$sum/$n}]
    if {$n >= 2} {
	set err [expr {sqrt(($sumSq - $sum*$sum/$n)/($n-1)/$n)}] 
    } else {
	set err 0.0
    }

    # Write the result.
    puts $out "$tim $mean $err"
    
    close $out
}
