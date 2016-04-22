# Shift a column in the data.
# Use with: tclsh addColumn.tcl 0 shift input_file
# Author: Jeff Comer <jcomer2@illinois.edu>

if {$argc != 1} {
	puts "This script requires one argument: "
	puts "  input file"
	puts ""
	exit
}

# Input:
set inData [lindex $argv 0]
# Output:
set outData ${inData}.stress

# Determine the number of sections.
set in [open $inData r]
set first 0
set sections 0
foreach line [split [read $in] \n] {
	if {[string match "STRESS*" $line]} {
		if {$first == 1} {break}
		set first 1
	} elseif {![string match "*AVGFx*" $line]} {
		set tok [concat $line]
		puts "Force entry: [lindex $tok 0]"
		set secArray([lindex $tok 0]) $sections
		incr sections
	}
}
puts "Reading $sections force entries."
close $in

# Open `sections' output files.
set out {}
for {set i 0} {$i < $sections} {incr i} {
	lappend out [open ${outData}${i} w]
}


# Write the data into the appropriate files.
set step 0
set in [open $inData r]
foreach line [split [read $in] \n] {	
	if {[string match "STRESS*" $line]} {
		set tok [concat $line]
		set step [lindex $tok 3]
	} elseif {![string match "*AVGFx*" $line]} {
		set tok [concat $line]
		if {[llength $tok] < 2} {continue}
		
		# This is numeric data.
		set sec $secArray([lindex $tok 0])
		
		# Select the output file.
		set o [lindex $out $sec]
		
		puts -nonewline $o $step
		for {set i 1} {$i < [llength $tok]} {incr i} {
			puts -nonewline $o " [lindex $tok $i]"
		}
		puts $o ""
	}
}
close $in

set i 0
foreach o $out {
	puts "The results were successfully written to `${outData}${i}'."
	incr i
	close $o
}






