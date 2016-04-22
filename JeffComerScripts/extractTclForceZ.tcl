# Compute block averages and standard deviations.
# Use with: tclsh blockAvgSe.tcl window_size input_size
# Author: Jeff Comer <jcomer2@uiuc.edu>

if {$argc < 2} {
    puts "Usage: tclsh blockAvgSe.tcl timestep inputFile0 \[inputFile1\]..."
    exit
}

set timestep [lindex $argv 0]
# Input:
set inFileList [lrange $argv 1 end]
# Output:
set outSuffix force

if {![string is double $timestep]} {
    puts "ERROR: The timestep size must be a real number."
    exit
}

set outFile [lindex $inFileList end].$outSuffix
set out [open $outFile w]

set n 0
set tLast 0.0
foreach inFile $inFileList {
    # Open the files.
    set in [open $inFile r]
    
    puts "Reading $inFile..."
    
    while {[gets $in line] >= 0} {
	if {[string length $line] <= 1} {continue}
	
	if {[string match "#*" $line]} {
	    puts $out $line
	} elseif {[string match "*TCLDNAFORCE*" $line]} {
	    set item [concat $line]
	    set t [expr $timestep*[lindex $item 2] + $tLast]
	    set f [lindex $item 6]

	    puts $out "$t $f"
	    incr n
	}

    }
    set tLast $t
    close $in
}

close $out
puts "$n forces were written to `${outFile}'."
