#!/usr/bin/tclsh
# Compute -grad(y) using a center difference.
# Author: Jeff Comer <jcomer2@illinois.edu>

if {$argc != 2} {
    puts "Usage: tclsh $argv0 inputDataFile xZeroLevel"
    exit
}
source $env(HOME)/scripts/useful.tcl

set inFile [lindex $argv 0]
set x0 [lindex $argv 1]


set data [readData $inFile]

set nearDist [expr {abs([lindex $data 0 0]-$x0)}]
set nearPot [expr {0.5*([lindex $data 0 1]+[lindex $data 1 1])/([lindex $data 0 0]-[lindex $data 1 0])}]

# Trapezoid integration.
set n [llength $data]
set sum 0.0
set potX {}
set potV {}
for {set i 1} {$i < $n} {incr i} {
    set p0 [lindex $data [expr {$i-1}]]
    set p1 [lindex $data $i]
    
    set dx [expr {[lindex $p1 0]-[lindex $p0 0]}]
    set mx [expr {0.5*([lindex $p1 0]+[lindex $p0 0])}]
    set my [expr {0.5*([lindex $p1 1]+[lindex $p0 1])}]

    set sum [expr {$sum + $dx*$my}]
    lappend potX $mx
    lappend potV $sum

    set dist [expr {abs($mx-$x0)}]
    if {$dist < $nearDist} {
	set nearDist $dist
	set nearPot $sum
    }
}

# Set the zero and print the output.
foreach x $potX v $potV {
    puts "$x [expr {$nearPot-$v}]"
}

exit
