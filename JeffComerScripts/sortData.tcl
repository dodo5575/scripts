#!/usr/bin/tclsh
# Author: Jeff Comer <jcomer2@illinois.edu>

if {$argc < 2} {
    puts "Usage: tclsh blockAvgSe.tcl column inputFile0 \[inputFile1\]..."
    exit
}

set column [lindex $argv 0]
# Input:
set inFileList [lrange $argv 1 end]
# Output:
set outSuffix $column.sorted

puts "Sorting rows by column ${column}..."
if {![string is integer $column]} {
    puts "ERROR: The column must be an integer."
    exit
}

foreach inFile $inFileList {
    # Open the files.
    set in [open $inFile r]
    puts "Sorting $inFile..."
    
    set data {}
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
		if {$nCol <= $column} {
		    puts "ERROR: Column $column does not exist."
		    exit
		}
	    } elseif {[expr $nCol1 != $nCol]} {
		puts "Warning: Expecting $nCol data columns but found $nCol1!"
		continue
	    }

	    lappend data $item
	}
    }
    close $in
    
    set outFile $inFile.$outSuffix
    #set out [open $outFile w]
    foreach row [lsort -real -index $column $data] {
	#puts $out $row
	puts $row
    }
    
    #close $out
    puts "`${outFile}' has been sorted."
}



