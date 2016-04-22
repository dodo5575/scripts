#!/usr/bin/tclsh
# Compute -grad(y) using a center difference.
# Author: Jeff Comer <jcomer2@illinois.edu>

if {$argc != 1} {
    puts "Usage: tclsh $argv0 inputDataFile"
    exit
}
source $env(HOME)/scripts/useful.tcl

set data [readData [lindex $argv 0]]

set n [llength $data]
for {set i 1} {$i < $n} {incr i} {
    set p0 [lindex $data [expr {$i-1}]]
    set p1 [lindex $data $i]

    set x [expr {0.5*([lindex $p0 0] + [lindex $p1 0])}]
    set dx [expr {[lindex $p1 0]-[lindex $p0 0]}]
    set dy [expr {[lindex $p1 1]-[lindex $p0 1]}]
    set force [expr {-$dy/$dx}]
    puts "$x $force"
}

exit
