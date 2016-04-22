# Shift a column in the data.
# Use with: tclsh addColumn.tcl 0 shift input_file
# Author: Jeff Comer <jcomer2@illinois.edu>

if {$argc != 3} {
	puts "This script requires three arguments: "
	puts "  column (zero based)"
	puts "  shift"
	puts "  input file"
	puts ""
	exit
}

set column [lindex $argv 0]
set shift [expr 1.0*[lindex $argv 1]]
# Input:
set inData [lindex $argv 2]
# Output:
set outData ${inData}.shifted

# Write the header.
set out [open $outData w]
puts "Shifting column ${column} by ${shift}..."

# Get the data.
set data {}
set in [open $inData r]
set nCol -1
foreach line [split [read $in] \n] {
	if {[string match "#*" $line]} {
		puts $out $line
	} else {
		set item [split [string trim $line]]
		set nCol1 [llength $item]
		
		if {$nCol < 0} {
			set nCol $nCol1
			puts "Reading $nCol columns."
		        lappend data $item
		} elseif {[expr $nCol1 != $nCol]} {
		        if {$nCol1 != 0} {
			    puts "Warning: Expecting $nCol data columns but found $nCol1!"
		        }
		} else {
		        lappend data $item
		}
	}
}
close $in

# Process.
set result {}
foreach d $data {
	set row {}
	for {set i 0} {$i < $nCol} {incr i} {
		set x [lindex $d $i]
		if {$i == $column} {
			set x [expr $x+$shift]
		}
		lappend row $x
	}
	lappend result $row
}

foreach d $result {
	foreach x $d {
		puts -nonewline $out "$x "
	}
	puts $out ""
}
puts "The results were successfully written to `${outData}'."
close $out





