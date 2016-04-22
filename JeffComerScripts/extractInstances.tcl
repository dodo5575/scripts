#! /usr/bin/tclsh
# Author: Jeff Comer <jcomer2@illinois.edu>

set atomSet {{11916 11917 11918 11919 11920 11921 11922} {11911 11912 11913 11914 11915}}
set nameSet {pot chl}
# Input:
# Output:
set outName inst_at

if {$argc < 1} {
    puts "Usage: tclsh extractInstances.tcl inputFile0 \[inputFile1\]..."
    exit
}
set inFileList $argv

set out {}
foreach name $nameSet {
    lappend out [open ${outName}_${name}.dat w]
}

proc getType {atom} {
    global atomSet
    
    set ind 0
    foreach s $atomSet {
	if {[lsearch $s $atom] >= 0} {
	    return $ind
	}
	incr ind
    }
    puts "Error! Could not find index $atom."
    exit
}

foreach inFile $inFileList {
    # Open the files.
    set in [open $inFile r]
    puts "Read $inFile."
    set count 0

    set lineList {}
    while {[gets $in line] >= 0} {
	if {[string length $line] <= 1} { continue }
	if {[string match "#*" $line]} { continue }
	
	set item [concat $line]
	if {[llength $item] < 9} { 
	    puts "Only [llength $item] of 9 columns found."
	    continue 
	}

	lappend lineList $line
    }
    close $in

    if {[llength $lineList] <= 0} { continue }
    set lineSort [lsort -integer -index 0 $lineList]
    set inst0 [lindex $lineSort 0 1]
    set atom0 [lindex $lineSort 0 0]
    set type [getType $atom0]

    foreach l $lineSort {
	set atom [lindex $l 0]
	set inst [lindex $l 1]

	if {$atom != $atom0} {
	    puts [lindex $out $type] "END"
	    set atom0 $atom
	    set inst0 $inst
	    set type [getType $atom]
	    incr count
	} elseif {$inst != $inst0} {
	    puts [lindex $out $type] "END"
	    set inst0 $inst
	    incr count
	}

	puts [lindex $out $type] [lrange $l 2 end]
    }

    puts "Found $count instances."
    if {$count > 0} {
	puts [lindex $out $type] "END"
    }
}

foreach o $out {
    close $o
}
puts "Complete."
exit
