# Change a type
# Author: Jeff Comer <jcomer2@illinois.edu>

if {$argc != 3} {
	puts "This script requires three arguments: "
	puts "  input type"
	puts "  output type"
	puts "  input file"
	puts ""
	exit
}

set inType [lindex $argv 0]
set outType [lindex $argv 1]
# Input:
set inData [lindex $argv 2]
# Output:
set outData ${inData}.type

if {[file exists $outData]} {
    puts "Error: File $outData exists!"
    exit
}

set in [open $inData r]
set out [open $outData w]
set count 0
foreach line [split [read $in] \n] {
    if {[string match "!*" $line]} {continue}

    set n [regexp -all $inType $line]
    if {$n >= 0} {
	# Replace each instance of inType with outType one at a time.
	# Hence, we assume not interaction between outType atoms.
	set start 0
	for {set i 0} {$i < $n} {incr i} {
	    regexp -start $start -indices $inType $line ind
	    puts $out [regsub -start $start $inType $line $outType]
	    set start [expr [lindex $ind 1]+1]
	}

	incr count
    }
}

close $out
close $in

puts "Extracted $count lines."
puts "Wrote `$outData'."



