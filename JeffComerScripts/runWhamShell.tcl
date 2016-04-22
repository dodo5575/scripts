# runWhamPhant.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

if {$argc > 0} {
    set name [lindex $argv 0]
} else {
    set name chl-chl
}

set indexFile three_${name}_ions.txt
set outFile shell_reprise_${name}
set binN 50
set tol 0.001
set dir .
set cut0 0.1
set cut1 1e100
#set levelZ0 9
#set levelZ1 11
set kbtToKcalMol 0.5862292

source 1Dwham.tcl

proc readIndex {fileName} {
    set in [open $fileName r]
    set index {}

    while {[gets $in line] >= 0} {
	if {[string length $line] <= 1} { continue }
	if {[string equal [string index $line 0] "\#"] } { continue }
	
	lappend index [concat $line]
    }
    
    close $in
    return $index
}

proc readZ {fileName cut0 cut1} {
    set in [open $fileName r]
    set z {}
    
    while {[gets $in line] >= 0} {
	if {[string length $line] <= 1} { continue }
	if {[string equal [string index $line 0] "\#"] } { continue }
	
	set tok [concat $line]
	set t [lindex $tok 0]
	if {[llength $tok] >= 2 && $t >= $cut0 && $t < $cut1} {
	    lappend z [lindex $tok 3]
	}
    }
    
    close $in
    return $z
}

#set hist [readHistogram $histFile]
set index [readIndex $indexFile]
set z {}
set cen {}
set spring {}

# Gather the data.
foreach rec $index {
    set dataFile "$dir/[lindex $rec 0]"
    lappend cen [lindex $rec 3]
    lappend spring [expr [lindex $rec 6]/$kbtToKcalMol]

    set data [readZ $dataFile $cut0 $cut1]
    lappend z $data

    puts "Read [llength $data] points from $dataFile for the window at z=[lindex $rec 3] with k=[lindex $rec 6]."
}

set levelZ0 [lindex $cen end-3]
set levelZ1 [lindex $cen end-1]

set ret [wham z $binN $cen $spring $outFile $levelZ0 $levelZ1 $tol shell]

exit
