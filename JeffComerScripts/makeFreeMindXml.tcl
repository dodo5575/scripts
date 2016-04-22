#!/usr/bin/tclsh
# Author: Jeff Comer <jcomer2@illinois.edu>

if {$argc != 1} {
    puts "Usage: tclsh $argv0 inputFile0"
    exit
}

set inFile [lindex $argv 0]
# Output:
set outName $inFile.mm

# Open the files.
set in [open $inFile r]
set out [open $outName w]
puts "Converting $inFile into FreeMind xml format."

set nCol -1
while {[gets $in line] >= 0} {
    if {[string length $line] <= 1} { continue }

    if {[string match "#*" $line]} { continue }
    
    if {[string match "-*" $line]} {
	set name [string range $line 1 end]
	puts $out "<node COLOR='\#000000' POSITION='right' STYLE='fork' TEXT='$name'>"
    } elseif {[string match "/*" $line]} {
	set name [string range $line 1 end]
	puts $out "</node>"
    } else {
	set name $line
	puts $out "<node COLOR='\#000000' POSITION='right' STYLE='fork' TEXT='$name'/>"
    }
}

close $out
