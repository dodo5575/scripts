# Compute block averages and standard deviations.
# Use with: tclsh blockAvgSe.tcl window_size input_size
# Author: Jeff Comer <jcomer2@illinois.edu>

if {$argc < 2} {
    puts "Usage: tclsh blockAvgSe.tcl dt inputFile0 \[inputFile1\]..."
    exit
}

set dt [lindex $argv 0]
# Input:
set inFileList [lrange $argv 1 end]
# Output:
set outSuffix traj

puts "The timestep is $dt."

foreach inFile $inFileList {
    # Open the files.
    set in [open $inFile r]
    set outFile $inFile.$outSuffix
    set out [open $outFile w]
    puts "Adding time to $inFile..."
        
    # Initialize the sums.
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
	    } elseif {[expr $nCol1 != $nCol]} {
		puts "Warning: Expecting $nCol data columns but found $nCol1!"
		continue
	    }
	    
	    set x [expr $n*$dt]
	    set y [lindex $item 0]
	    puts $out "$x $y"
	    incr n
	}
    }
 
    close $out
    close $in
    puts "The two columns were successfully written to `${outFile}'."
}



