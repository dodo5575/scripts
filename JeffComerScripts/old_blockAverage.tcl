# Make block averages.
# Use with: tclsh blockAverage.tcl window_size input_size
# Author: Jeff Comer <jcomer2@illinois.edu>

if {$argc != 2} {
	puts "This script requires two arguments: "
	puts "  window size"
	puts "  input file"
	puts ""
	exit
}

set window [expr 1.0*[lindex $argv 0]]
# Input:
set inData [lindex $argv 1]
# Output:
set outData ${inData}.block

# Write the header.
set out [open $outData w]
#puts $out "\# Block averages with window $window"
puts "Forming block averages with a window of size ${window}..."

# Get the data.
set data {}
set in [open $inData r]
set nCol -1
foreach line [split [read $in] \n] {
	if {[string length $line] <= 1} {break}
	
	if {[string match "#*" $line]} {
		puts $out $line
	} else {
		set item [split [string trim $line]]
		set nCol1 [llength $item]
		
		if {$nCol < 0} {
			set nCol $nCol1
			puts "Reading $nCol columns."		
		} elseif {[expr $nCol1 != $nCol]} {
			puts "ERROR: Expecting $nCol data columns but found $nCol1!"
			exit
		}
		lappend data $item
	}
}
close $in

# Create the zero row.
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
	puts "Making a final data point with leftover data..."
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
puts "The block averages were successfully written to `${outData}'."
close $out





