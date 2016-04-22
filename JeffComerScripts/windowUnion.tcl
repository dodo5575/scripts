#!/usr/bin/tclsh
# Author: Jeff Comer <jcomer2@uiuc.edu>

if {$argc != 2} {
    puts "Usage: tclsh $argv0 masterWindows subWindows"
    exit
}

source $env(HOME)/scripts/useful.tcl

set data0 [readData [lindex $argv 0]]
set data1 [readData [lindex $argv 1]]

set pos1List {}
foreach d1 $data1 {
    lappend pos1List [lrange $d1 1 3]
}

foreach d0 $data0 {
    set pos [lrange $d0 1 3]
    if {[lsearch $pos1List $pos] >= 0} {
	puts $d0
    }
}


