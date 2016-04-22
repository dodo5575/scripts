#! /usr/bin/tclsh
# Author: Jeff Comer <jcomer2@illinois.edu>

if {$argc != 1} {
    puts "Usage: tclsh dispBinary inFile"
    exit
}

set chunk 36
set inFile [lindex $argv 0]
set in [open $inFile r]
fconfigure $in -translation binary

set inBinData [read $in]
close $in
#puts "length [string length $inBinData]"

set good 1
while {[string length $inBinData] >= $chunk} {
    #puts $inData
    set count [binary scan $inBinData f3f3f3 pos0 pos1 force]
    puts "$pos0 $pos1 $force"

    set inBinData [string range $inBinData $chunk end]
}
