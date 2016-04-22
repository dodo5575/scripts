# Shift a column in the data.
# Use with: tclsh addColumn.tcl 0 shift input_file
# Author: Jeff Comer <jcomer2@illinois.edu>

if {$argc != 1} {
	puts "This script requires one arguments: "
	puts "  input file"
	puts ""
	exit
}

# Input:
set inData [lindex $argv 0]
# Output:
set outData ${inData}.nc

# Write the header.
puts "Removing comments from $inData"
set out [open $outData w]

# Get the data.
set data {}
set in [open $inData r]
set nCol -1
foreach line [split [read $in] \n] {
	if {![string match "#*" $line]} {
		puts $out $line
	}
}
close $in
puts "The results were successfully written to `${outData}'."
close $out





