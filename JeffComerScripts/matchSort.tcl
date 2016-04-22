# Author: Jeff Comer <jcomer2@illinois.edu>

# Input:
set inFile relative_helix3_100cM.txt
set matchFile relative_triplet_100cM.txt
# Output:
set sortFile relative_helix3_100cM_sort.txt

# Just read a space delimited data file.
proc readData {fileName} {
    set in [open $fileName r]
    
    set r {}
    while {[gets $in line] >= 0} {
	if {[string match "#*" $line]} {continue}
	if {[string length $line] < 2} {continue}

	set tok [concat $line]
	set t [lindex $tok 0]
	lappend r $tok
    }

    close $in
    return $r
}

set matchData [readData $matchFile]
set inData [readData $inFile]
puts "Read [llength $matchData] items from `$matchFile'."
puts "Read [llength $inData] items from `$inFile'."

set out [open $sortFile w]
set count 0
foreach m $matchData {
    set key [lindex $m 0]

    # Find this key in the input data and write it.
    foreach i $inData {
	if {[string equal $key [lindex $i 0]]} {
	    puts $out $i
	    incr count
	    break
	}
    }
}
close $out

puts "Matched $count keys."
