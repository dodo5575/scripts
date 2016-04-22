# Shift a column.
# Author: Jeff Comer <jcomer2@illinois.edu>

if {$argc < 2} {
    puts "Usage: tclsh removeBaseline.tcl intercept slope inputFile0 \[inputFile1\]..."
    exit
}

set intercept [lindex $argv 0]
set slope [lindex $argv 1]
# Input:
set inFileList [lrange $argv 2 end]
# Output:
set outSuffix shifted

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
	    if {$nCol < 2} {
		puts "ERROR: Two columns required."
		exit
	    }
	} elseif {[expr $nCol1 != $nCol]} {
	    puts "Warning: Expecting $nCol data columns but found $nCol1!"
	    continue
	}

	set x [lindex $item 0]
	set y [lindex $item 1]
	set y1 [expr $y - ($intercept + $slope*$x)]
	puts $out "$x $y1"
    }
    
    close $out
    close $in
    puts "Results were successfully written to `${outFile}'."
}



