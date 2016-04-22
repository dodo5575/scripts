# Shift a column in the data.
# Use with: tclsh addColumn.tcl 0 shift input_file
# Author: Jeff Comer <jcomer2@illinois.edu>

if {$argc != 3} {
	puts "This script requires three arguments: "
	puts "  column x (zero based)"
	puts "  column y (zero based)"
	puts "  input file"
	puts ""
	exit
}

set columnX [lindex $argv 0]
set columnY [lindex $argv 1]
# Input:
set inData [lindex $argv 2]
# Output:
set outData ${inData}.${columnX}_${columnY}.col

proc tokenList {str} {
	set ret {}
	set item [split $str " "]
	
	foreach i $item {
		if {[llength $i] != 0} {
			lappend ret $i
		}
	}
	
	return $ret
}

# Write the header.
set out [open $outData w]
puts "Extracting columns ${columnX} and ${columnY}"

# Get the data.
set in [open $inData r]
set nCol -1
foreach line [split [read $in] \n] {
	if {[string match "#*" $line]} {
		# Write any comments.
		puts $out $line
	} else {
		set item [tokenList $line]
		set nCol1 [llength $item]
		
		if {$nCol1 == 0} {
			continue
		}
	
		if {$nCol1 != $nCol && $nCol > 0} {
			puts "Warning: Expecting $nCol data columns but found $nCol1!"
			continue
		}
		
		if {$nCol < 0} {
			set nCol $nCol1
			puts "Reading $nCol columns."    
		}
		
		# Write the desired columns.
		puts $out "[lindex $item $columnX] [lindex $item $columnY]"
	}
}
close $in
puts "The results were successfully written to `${outData}'."
close $out





