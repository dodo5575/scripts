# Shift a column in the data.
# Use with: tclsh addColumn.tcl shift input_file
# Author: Jeff Comer <jcomer2@illinois.edu>

if {$argc != 2} {
    puts "This script requires two arguments: "
    puts "  sigma"
    puts "  input file"
    puts ""
    exit
}

set sigma [expr {1.0*[lindex $argv 0]}]
set sigmaSq [expr {$sigma*$sigma}]
# Input:
set inData [lindex $argv 1]
# Output:
set outData ${inData}.gauss${sigma}

# Write the header.
set out [open $outData w]
puts "Gaussian filter"

# Get the data.
set tData {}
set xData {}
set in [open $inData r]
set nCol -1
foreach line [split [read $in] \n] {
    if {[string length $line] <= 1} { break }
    
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
	lappend tData [expr {double([lindex $item 0])}]
	lappend xData [expr {double([lindex $item 1])}]
    }
}
close $in

# Process.
set n [llength $tData]
set dt [expr {([lindex $tData end] - [lindex $tData 0])/($n-1)}]
puts "Average data interval: $dt"
set window [expr {int(8.0*$sigma/$dt)}]
set windowHalf [expr {$window/2}]
puts "Using a sampling window of eight standard deviations: $window"

set result {}
set index 0
foreach t $tData {
    
    set sumWeights 0.0
    set sum 0.0
    
    set start [expr {$index - $windowHalf}]
    set stop [expr {$index + $windowHalf}]
    if {$start < 0} {set start 0}
    if {$stop >= $n} {set stop [expr {$n-1}]}
    
    for {set i $start} {$i < $stop} {incr i} {
	set ti [lindex $tData $i]
	set xi [lindex $xData $i]
	
	set wi [expr {exp(-($t-$ti)*($t-$ti)/(2.0*$sigmaSq))}]
	set sumWeights [expr {$sumWeights + $wi}]
	set sum [expr {$sum + $wi*$xi}]
    }
    set weightedMean [expr {$sum/$sumWeights}]
    
    lappend result [list $t $weightedMean]
    incr index
}

foreach r $result {
    puts $out "[lindex $r 0] [lindex $r 1]" 
}
puts "The results were successfully written to `${outData}'."
close $out





