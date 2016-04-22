# Make block averages.
# Use with: tclsh blockAverage.tcl window_size input_size
# Author: Jeff Comer <jcomer2@illinois.edu>

if {$argc != 2} {
	puts "This script requires two arguments: "
	puts "  stride"
	puts "  input file"
	puts ""
	exit
}

set stride [lindex $argv 0]
# Input:
set inData [lindex $argv 1]
# Output:
set outData ${inData}.skip${stride}

set in [open $inData r]
set out [open $outData w]
set count 0
foreach line [split [read $in] \n] {
	if {![expr $count%$stride]} {
		puts $out $line
	}
	incr count
}

close $out
close $in




