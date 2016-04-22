# Author: Jeff Comer <jcomer2@illinois.edu>

if {$argc < 2} {
    puts "Usage: tclsh $argv0 potentialFile potentialMax"
    exit
}
set potentialFile [lindex $argv 0]
set potentialMax [lindex $argv 1]

# Just read a space delimited data file.
proc readData {fileName} {
    set in [open $fileName r]
    
    set r {}
    while {[gets $in line] >= 0} {
	if {[string match "#*" $line]} {continue}
	if {[string length $line] < 2} {continue}

	set tok [concat $line]
	lappend r $tok
    }

    close $in
    return $r
}

set data [readData $potentialFile]
set dx [expr {[lindex $data 1 0]-[lindex $data 0 0]}]
set x0 [lindex $data 0 0]
set v0 [lindex $data 0 1]

set m [expr {($potentialMax-$v0)/(0.0-$x0)}]
set n [expr {int(floor($x0/$dx))}]

for {set i $n} {$i > 0} {incr i -1} {
    set x [expr {$x0 - $i*$dx}]
    set v [expr {$m * ($x - $x0) + $v0}]
    puts [list $x $v]
}

foreach d $data {
    puts $d
}



