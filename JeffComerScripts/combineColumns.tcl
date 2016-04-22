# Shift a column in the data.
# Use with: tclsh combineColumns inFile0 inFile1 outFile
# Author: Jeff Comer <jcomer2@illinois.edu>

if {$argc != 3} {
	puts "This script requires three arguments: "
	puts "  inFile0"
	puts "  inFile1"
	puts "  outFile"
	puts ""
	exit
}

set col0 1
set col1 1
# Input:
set inFile0 [lindex $argv 0]
set inFile1 [lindex $argv 1]
# Output:
set outFile [lindex $argv 2]

# Write the header.
puts "Combining column $col0 from $inFile0 and column $col1 from $inFile1 together in $outFile."

# Get the data.
set data0 {}
set in [open $inFile0 r]
puts "Reading $inFile0..."
set nCol -1
foreach line [split [read $in] \n] {
    if {[string length $line] <= 1} {break}
    
    if {[string match "#*" $line]} {
	puts $out $line
    } else {
	set item [concat $line]
	set nCol1 [llength $item]
	if {$nCol < 0} {
	    set nCol $nCol1
	    puts "Reading $nCol columns."
	} elseif {[expr $nCol1 != $nCol]} {
	    puts "ERROR: Expecting $nCol data columns but found $nCol1!"
	    exit
	}
	lappend data0 [lindex $item $col0]
    }
}
close $in

# Get the data.
set data1 {}
set in [open $inFile1 r]
puts "Reading $inFile1..."
set nCol -1
foreach line [split [read $in] \n] {
    if {[string length $line] <= 1} {break}
    
    if {[string match "#*" $line]} {
	puts $out $line
    } else {
	set item [concat $line]
	set nCol1 [llength $item]
	if {$nCol < 0} {
	    set nCol $nCol1
	    puts "Reading $nCol columns."
	} elseif {[expr $nCol1 != $nCol]} {
	    puts "ERROR: Expecting $nCol data columns but found $nCol1!"
	    exit
	}
	lappend data1 [lindex $item $col1]
    }
}
close $in

# Output.
set out [open $outFile w]
foreach d0 $data0 d1 $data1 {
    puts $out "$d0 $d1"
}

puts "The results were successfully written to `${outFile}'."
close $out



