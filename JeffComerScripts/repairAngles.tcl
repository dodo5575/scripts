# Compute block averages and standard deviations.
# Use with: tclsh blockAvgSe.tcl window_size input_size
# Author: Jeff Comer <jcomer2@illinois.edu>

if {$argc < 2} {
    puts "Usage: tclsh repairAngles.tcl window inputFile0 [inputFile1]..."
    exit
}

set window [lindex $argv 0]
set theta 360.0
set halfTheta [expr 0.5*$theta]
# Input:
set inFileList [lrange $argv 1 end]
# Output:
set outSuffix ang

set w [expr abs(int($window/2))]
puts "Repairing angles with window of [expr 2*$w]..."

foreach inFile $inFileList {
    # Get the data.
    set tData {}
    set xData {}
    set xs 0.0
    set in [open $inFile r]
    puts "Reading $inFile..."

    set nCol -1
    foreach line [split [read $in] \n] {
	if {[string length $line] <= 1} {break}
	
	if {[string match "#*" $line]} {
	    puts $out $line
	} else {
	    set item [concat $line]
	    set nCol1 [llength $item]
	    if {$nCol < 0} {
		set nCol $nCol1
		puts "Reading $nCol columns."
	    } elseif {[expr $nCol1 != $nCol]} {
		puts "ERROR: Expecting $nCol data columns but found $nCol1!"
		exit
	    }
	    lappend tData [lindex $item 0]
	    lappend xData [lindex $item 1]
	    set xs [expr $xs + [lindex $item 1]]
	}
    }
    close $in
    set xm [expr $xs/[llength $xData]]
    puts "Mean: $xm"

    set outFile $inFile.$outSuffix
    set out [open $outFile w]
    # Shift by theta if the deviation of x from the median of its block
    # is larger than theta/2.
    set n [llength $xData]
    set i 0
    foreach t $tData x $xData {
	set j0 [expr $i-$w]
	set j1 [expr $i+$w]
	if {$j0 < 0} {set j0 0}
	if {$j1 >= $n} {set j1 [expr $n-1]}

	# Find the median around this point.
	set region [lsort -real [lrange $xData $j0 $j1]]
	set mid [expr ($j1-$j0)/2]
	set xMedian [lindex $region $mid]
		
	if {$x - $xMedian > $halfTheta} {set x [expr $x - $theta]}
	if {$x - $xMedian < -$halfTheta} {set x [expr $x + $theta]}

	puts $out "$t $xMedian"
	incr i
    }
    close $out
}



