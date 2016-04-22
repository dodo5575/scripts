# Rename a set of segments.
# Use with: tclsh templateGraph.tcl templateFile [inputFile0] [inputFile1] ... outputFile
# Author: Jeff Comer <jcomer2@illinois.edu>

if {$argc != 3} {
    puts "This script requires at least two arguments: "
    puts "  templateFile"
    puts "  outputFile"
    puts "usage:  tclsh templateGraph.tcl templateFile inputFile0 outputFile"
    
    exit
}

set lastArg [expr $argc - 1]
set lastInp [expr $argc - 2]
# Input:
set templateFile [lindex $argv 0]
set inputFile [lrange $argv 1 $lastInp]
# Output:
set outputFile [lindex $argv $lastArg]

puts ""
puts "labelGraph:"
puts "Using $templateFile as a template"
puts "Using $inputFile as input data sets"
puts "Inserting [llength $inputFile] data sets into $outputFile"

if {[file exists $outputFile]} {
    puts stderr "ERROR: Output file $outputFile exists!"
    exit
}

# Insert a data file into an output stream.
proc insertData {out file} {
    set in [open $file r]
    set data [split [read $in] \n]
    close $in

    set len [expr {[llength $data]/2}]

    puts $out "@    xaxis  tick spec type both"
    puts $out "@    xaxis  tick spec $len"
    
    foreach item $data {
	if {[string length $item] < 1} { continue }

	puts $out $item
    }
}

set in [open $templateFile r]
set out [open $outputFile w]
set count 0
set readingData 0
foreach line [split [read $in] \n] {
    if {$readingData} {
	puts $out $line
	set readingData 0
    } elseif {[string match "@    xaxis  tick spec*" $line] && $count == 0} {
	# Begin reading a data set.
	set readingData 1
	# Insert the data file.
	if {$count >= [llength $inputFile]} {
	    puts "Warning: No input file to insert. Data set $count will be blank."
	    incr count
	} else {
	    set current [lindex $inputFile $count]
	    puts "Data set $count: $current"
	    insertData $out $current
	    incr count
	}
    } else {
	# Just write the line.
	puts $out $line
    }

	  }
close $out
close $in

puts ""
puts "Found $count insertion points."
puts "Had [llength $inputFile] sets to insert."
exit



