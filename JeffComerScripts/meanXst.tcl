# Author: Jeff Comer <jcomer2@illinois.edu>

if {$argc < 1} {
    puts "Usage: tclsh readXst inputFile0 \[inputFile1\]..."
    exit
}

# Input:
set inFileList [lrange $argv 0 end]
# Output:
set outSuffix dat


set dimList {cz}
foreach inFile $inFileList {
    # Open the files.
    set in [open $inFile r]
        
    foreach dim $dimList {
	set out($dim) [open $inFile.$dim.$outSuffix w]
	set sum($dim) 0
	set sumSq($dim) 0
    }
    set n 0

    while {[gets $in line] >= 0} {
	if {[string length $line] <= 1} {break}
	
	if {[string match "#*" $line]} {continue}
	set tok [concat $line]
	set x(ax) [lindex $tok 1]
	set x(ay) [lindex $tok 2]
	set x(bx) [lindex $tok 4]
	set x(by) [lindex $tok 5]
	set x(cz) [lindex $tok 9]

	foreach dim $dimList {
	    puts $out($dim) "[lindex $tok 0] $x($dim)"
	    set sum($dim) [expr $sum($dim) + $x($dim)]
	    set sumSq($dim) [expr $sumSq($dim) + $x($dim)*$x($dim)]
	}
	
	incr n
    }
    
    puts "\nRead $inFile"
    foreach dim $dimList {
	set mean($dim) [expr $sum($dim)/$n]
	set err($dim) [expr sqrt(($sumSq($dim) - $sum($dim)*$sum($dim)/$n)/($n-1)/$n)]
	puts "mean $dim: $mean($dim) +/- $err($dim)"
    }
    
    foreach dim $dimList {
	close $out($dim)
    }
    close $in
}

