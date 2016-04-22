# Remove lines matching a pattern.
# to use: tclsh removeLines.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set pattern "*2.0\**"

# Input:
set inFile diffuse_DS3_fast.txt
# Output:
set outFile diffuse_DS3_K.txt

set in [open $inFile r]
set out [open $outFile w]

foreach line [split [read $in] "\n"] {
	if {![string match $pattern $line]} {
		puts $out $line
	}
}

close $out
close $in




