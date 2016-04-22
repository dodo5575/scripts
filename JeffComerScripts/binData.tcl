# Compute block averages and standard deviations.
# Use with: tclsh blockAvgStd.tcl window_size input_size
# Author: Jeff Comer <jcomer2@illinois.edu>

if {$argc < 3} {
    puts "Usage: tclsh blockAvgStd.tcl windowSize inputFile0 \[inputFile1\]... outFile"
    exit
}

set window [lindex $argv 0]
# Input:
set inFileList [lrange $argv 1 end-1]
# Output:
set outFile [lindex $argv end]

puts "Forming block averages with a window of size ${window}..."
if {![string is double $window]} {
    puts "ERROR: The window size must be a number."
    exit
}


foreach inFile $inFileList {
    # Open the files.
    set in [open $inFile r]
    puts "Blocking $inFile..."
    set n 0
    
    set nCol -1
    while {[gets $in line] >= 0} {
	if {[string length $line] <= 1} {break}
	
	if {[string match "#*" $line]} {
	    puts $out $line
	} else {
	    set item [concat $line]
	    set nCol1 [llength $item]
	    if {$nCol < 0} {
		set nCol $nCol1
		puts "Reading $nCol columns."
		if {$nCol < 2} {
		    puts "ERROR: At least two columns are required."
		    exit
		}
		if {$nCol != 2} {puts "Warning: Ignoring [expr $nCol-2] columns."}
	    } elseif {[expr $nCol1 != $nCol]} {
		puts "Warning: Expecting $nCol data columns but found $nCol1!"
		continue
	    }
	    
	    # Bin.
	    incr n
   
	    set x [lindex $item 0]
	    set y [lindex $item 1]
	    set b [expr int(floor($x/$window))]
	    lappend bin(b${b}) [list $x $y]
	}
    }
    close $in
    puts "$n values were read."
}

set binAvg {}
puts "Forming the bin averages..."
foreach {binName binList} [array get bin] {
    set n [llength $binList]
    if {$n <= 1} {continue}
    
    set xSum 0
    set ySum 0
    set xxSum 0
    set yySum 0

    foreach pair $binList {
	foreach {x y} $pair {break}
	set xSum [expr $xSum + $x]
	set ySum [expr $ySum + $y]
	set xxSum [expr $xxSum + $x*$x]
	set yySum [expr $yySum + $y*$y]
    }

    set xAvg [expr $xSum/$n]
    set yAvg [expr $ySum/$n]
    set xErr [expr sqrt(($xxSum + $xSum*$xSum/$n)/($n-1)/$n)]
    set yErr [expr sqrt(($yySum + $ySum*$ySum/$n)/($n-1)/$n)]

    lappend binAvg [list $xAvg $yAvg $yErr]

    puts "Bin: $xAvg $n"
}

# Write the results.
set out [open $outFile w]
foreach b [lsort -index 0 -real $binAvg] {
    puts $out $b
}
close $out
puts "Bin averages were written to `${outFile}'."



