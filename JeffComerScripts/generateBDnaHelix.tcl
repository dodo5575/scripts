#!/usr/bin/tclsh
# Author: jcomer2@illinois.edu
if {$argc != 6} {
    puts "Usage: $argv0 sysGrid.dx bpGrid.dx bpNum bpRise bpPerTurn outputFile"
    exit
}

set sysGrid [lindex $argv 0]
set bpGrid [lindex $argv 1]
set bpNum [lindex $argv 2]
set bpRise [lindex $argv 3]
set bpPerTurn [lindex $argv 4]
set outFile [lindex $argv [expr $argc-1]]

if {![file exists $sysGrid]} {
    puts "ERROR: System grid `$sysGrid' does not exist."
    exit
}
if {![file exists $bpGrid]} {
    puts "ERROR: Basepair grid `$bpGrid' does not exist."
    exit
}


# Get the geometry.
set twoPi [expr {8.0*atan(1.0)}]
set theta [expr {-$twoPi/$bpPerTurn}]
set z0 [expr {-0.5*$bpNum*$bpRise}]

# Open the file and write the first line.
set out [open $outFile w]
puts $out "$sysGrid 1.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0"

for {set i 0} {$i < $bpNum} {incr i} {
    set ct [expr {cos($i*$theta)}]
    set st [expr {sin($i*$theta)}]
    set nst [expr {-$st}]

    set m "$ct $st 0.0 $nst $ct 0.0 0.0 0.0 1.0"
    set d "0.0 0.0 [expr {$z0 + $bpRise*$i}]"

    puts $out "$bpGrid $m $d"
}
close $out
