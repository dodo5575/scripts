# Compute block averages and standard deviations.
# Use with: tclsh blockAvgSe.tcl window_size input_size
# Author: Jeff Comer <jcomer2@illinois.edu>

if {$argc < 2} {
    puts "Usage: tclsh blockAvgSe.tcl windowSize inputFile0 [inputFile1]..."
    exit
}

set window [lindex $argv 0]
# Input:
set inFileList [lrange $argv 1 end]
# Output:
set outSuffix $window.block

puts "Forming block averages with a window of size ${window}..."
if {![string is integer $window]} {
    puts "ERROR: The window size must be an integer."
    exit
}

foreach inFile $inFileList {
    # Get the data.
    set tData {}
    set xData {}
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
	}
    }
    close $in

    # Form the block averages and standard errors.
    set tAvg {}
    set xAvg {}
    set xSe {}
    set n 0
    set tSum 0
    set xSum 0
    set xSumSq 0
    foreach t $tData x $xData {
	incr n

	set tSum [expr $tSum + $t]
	set xSum [expr $xSum + $x]
	set xSumSq [expr $xSumSq + $x*$x]

	if {$n == $window} {
	    set ta [expr $tSum/$n]
	    set xa [expr $xSum/$n]
	    set xxa [expr $xSumSq/$n]

	    lappend tAvg $ta
	    lappend xAvg $xa
	    lappend xSe [expr sqrt(($xxa - $xa*$xa)/$n)]

	    set n 0
	    set tSum 0
	    set xSum 0
	    set xSumSq 0
	}
    }

    # Make a final data point if we have a lot of data left.
    if {[expr $n > $window/2 + 1]} {
	puts "Making a final block average with the final $n data points..."
	set ta [expr $tSum/$n]
	set xa [expr $xSum/$n]
	set xxa [expr $xSumSq/$n]

	lappend tAvg $ta
	lappend xAvg $xa
	lappend xSe [expr sqrt(($xxa - $xa*$xa)/$n)]
    }

    # Write the block averages.
    set outFile $inFile.$outSuffix
    set out [open $outFile w]
    foreach t $tAvg x $xAvg s $xSe {
	puts $out "$t $x $s"
    }
    puts "The block averages and standard errors were successfully written to `${outFile}'."
    close $out
}



