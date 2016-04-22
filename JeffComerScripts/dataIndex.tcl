# This script generates the data mining index file used
# the Matlab script.
# It alerts you if any files fail to exist.
# Author: Jeff Comer <jcomer2@illinois.edu>

set name pores
# Input:
set inFile ${name}_index.txt
# Output:
set outFile ${name}_data.txt

set in [open $inFile r]
set out [open $outFile w]

foreach line [split [read $in] \n] {
    # Ignore comments.
    if {[string match "\#*" $line]} {continue}

    set tok [concat $line]
    if {[llength $tok] < 11} {continue}

    set class [lindex $tok 1]
    set name [lindex $tok 2]
    set range [lindex $tok 8]
    set offset [lindex $tok 10]
    set voltage  [string range [regexp -inline {\_[0-9.]*V} $name] 1 end-1]
    set fileName $name[lindex $range 0]-[lindex $range end]

    puts $out "$fileName $voltage $class $offset"
    puts $fileName
}

close $out
close $in



