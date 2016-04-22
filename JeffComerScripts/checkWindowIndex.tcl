#!/usr/bin/tclsh
# Author: Jeff Comer <jcomer2@uiuc.edu>

if {$argc != 2} {
    puts "Usage: tclsh $argv0 indexFile0 outFile"
    exit
}

source $env(HOME)/scripts/useful.tcl

set data [readData [lindex $argv 0]]
set outFile [lindex $argv 1]

set out [open $outFile w]

foreach d $data {
    set f [lindex $d 0]
    if {[file exists $f]} {
	puts $out $d
    } else {
	puts "File `$f' does not exist."
    }
}

close $out
