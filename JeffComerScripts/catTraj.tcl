#!/usr/bin/tclsh
# Concatenate time series files.
# Use with: tclsh catSeries.tcl inputFile0 inputFile1 inputFile2...
# Author: Jeff Comer <jcomer2@illinois.edu>

if {$argc < 1} {
	puts "Use with: tclsh catSeries.tcl inputFile0 inputFile1 inputFile2..."
	exit
}

set tim0 0.0
foreach inputFile $argv {
    set inp [open $inputFile r]
    while {[gets $inp line] >= 0} {
	if {[string match "\#*" $line]} { continue }

	set tok [concat $line]
	if {[llength $tok] < 2} {
	    puts $line
	    continue
	}
	set t [lindex $tok 1]
	puts "[lindex $tok 0] [expr {$t+$tim0}] [lrange $tok 2 end]"
    }
    close $inp
    # Set the new time offset.
    set tim0 [expr {$t+$tim0}]
}
