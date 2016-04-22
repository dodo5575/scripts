# Shift a column.
# Author: Jeff Comer <jcomer2@illinois.edu>

if {$argc < 2} {
    puts "Usage: tclsh addColumn.tcl column shift inputFile0 \[inputFile1\]..."
    exit
}

set column [lindex $argv 0]
set shift [lindex $argv 1]
# Input:
set inFileList [lrange $argv 2 end]
# Output:
set outSuffix shifted

puts "Shifting column $column..."
if {![string is integer $column]} {
    puts "ERROR: The column index (zero-based) must be an integer."
    exit
}

foreach inFile $inFileList {
    # Open the files.
    set in [open $inFile r]
    set outFile $inFile.$outSuffix
    set out [open $outFile w]
    puts "Shifting $inFile..."

    set nCol -1
    set last -1
    while {[gets $in line] >= 0} {
	if {[string match "#*" $line]} {
	    puts $out $line
	    continue
	}

	set item [concat $line]
	set nCol1 [llength $item]
	if {$nCol < 0} {
	    set nCol $nCol1
	    set last [expr $nCol-1]
	    puts "Reading $nCol columns."
	    if {$nCol <= $column} {
		puts "ERROR: Column $column selected, but data has only $nCol columns."
		exit
	    }
	} elseif {[expr $nCol1 != $nCol]} {
	    puts "Warning: Expecting $nCol data columns but found $nCol1!"
	    continue
	}

	for {set i 0} {$i < $last} {incr i} {
	    if {$i == $column} {
		puts -nonewline $out "[expr [lindex $item $i] + $shift] "
	    } else {
		puts -nonewline $out "[lindex $item $i] "
	    }
	}
	if {$i == $column} {
	    puts $out "[expr [lindex $item $i] + $shift] "
	} else {
	    puts $out "[lindex $item $i] "
	}
    }
    
    close $out
    close $in
    puts "Results were successfully written to `${outFile}'."
}



