# Make block averages.
# Use with: tclsh blockAverage.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set window 20.0
# Input:
set inData cond_conform_4V_good.csv
# Output:
set outData cond_conform_4V_block.csv

# Write the header.
set out [open $outData w]
puts $out "\# Block averages with window $window"
puts "Forming block averages with window ${window}..."

# Get the data.
set data {}
set in [open $inData r]
foreach line [split [read $in] \n] {
	if {[string length $line] <= 1} {break}
	
	if {[string match "#*" $line]} {
		puts $out $line
	} else {
		lappend data [split $line]
	}
}
close $in

# Create the zero row.
set nCol [llength [lindex $data 0]]
puts "Reading $nCol columns."
set rowZero {}
for {set i 0} {$i < $nCol} {incr i} {
	lappend rowZero 0.
}

# Form the block averages.
set result {}
set rowSum $rowZero
set n 0
foreach d $data {
	set rowSum0 $rowSum
	set rowSum {}
	foreach x $d xs $rowSum0 {
		lappend rowSum [expr $xs + $x]
	}
		
	incr n
	
	# Compute the block average if we've finished a window.
	if {$n == $window} {
		set rowAvg {}
		foreach xs $rowSum {
			lappend rowAvg [expr $xs/$window]
		}
		lappend result $rowAvg
				
		# Clear the sums.
		set rowSum $rowZero
		set n 0
	}
}

# Make a final data point if we have a lot of data left.
if {$n > $window/2} {
	puts "Making a final data point..."
	set rowAvg {}
	foreach xs $rowSum {
		lappend rowAvg [expr $xs/$n]
	}
	lappend result $rowAvg
}

# Write the blocks.
foreach d $result {
	foreach x $d {
		puts -nonewline $out "$x "
	}
	puts $out ""
}
puts "The block averages have been written successfully."
close $out





