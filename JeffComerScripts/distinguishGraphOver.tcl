# Author: Jeff Comer <jcomer2@illinois.edu>

set inFile triplet140cM_currents.txt
set outPrefix over_triplet140cM

proc complement {seq} {
    set s [string map {- + + - a t t a g c c g} $seq]
    return [string index $s 0][reverse [string range $s 1 3]]
}

proc reverse {s} {
    set n [string length $s]
    set ret ""

    for {set i 0} {$i < $n} {incr i} {
	set c [string index $s $i]
	set ret ${c}${ret}
    }
    return $ret
}

# Open the file.
set in [open $inFile r]

# Read the data.
set seqList {}
while {[gets $in line] >= 0} {
    if {[string length $line] <= 1} { continue }
    set tok [concat $line]

    if {[llength $tok] != 3} { continue }
    
    lappend seqList $tok
}
close $in

set overList {}
# Add the non-distinct partners.
foreach seq $seqList {
    foreach {s i di} $seq { break }
    lappend overList $seq
    lappend overList [list [complement $s] $i $di]
}

# Sort the currents.
set sortList [lsort -real -index 1 $overList]
set togetherList [lsort -index 0 $overList]

set n [llength $togetherList]
set combinedList {}
for {set i 0} {$i < $n} {incr i} {
    if {$i < $n-1 && [lindex $togetherList $i 0] == [lindex $togetherList [expr {$i+1}] 0]} {
	set name [lindex $togetherList $i 0]
	set i0 [lindex $togetherList $i 1]
	set i1 [lindex $togetherList [expr {$i+1}] 1]
	set e0 [lindex $togetherList $i 2]
	set e1 [lindex $togetherList [expr {$i+1}] 2]

	set mean [expr {0.5*($i0 + $i1)}]
	set err [expr {0.5*sqrt($e0*$e0 + $e1*$e1)}]
	lappend combinedList [list $name $mean $err]
	incr i
    } else {
	lappend combinedList [lindex $togetherList $i]
    }
}
set combinedList [lsort -real -index 1 $combinedList]

set uniList {}
for {set i 0} {$i < [llength $combinedList]} {incr i} {
    lappend uniList [lindex $combinedList $i]
    set n0 [lindex $combinedList $i 0]
    set n1 [lindex $combinedList [expr {$i+1}] 0]
    if {$i < $n-1 && [string equal $n0 [complement $n1]]} {
	incr i
    } 
}

# Write the results.
set out [open ${outPrefix}_index.dat w]
set outSort [open ${outPrefix}_sort.txt w]
set i 0
foreach sort $sortList {
    puts $out "$i [lindex $sort 1] [lindex $sort 2]"
    puts $outSort "$sort"
    incr i
}
close $outSort
close $out

set out [open ${outPrefix}_together.txt w]
foreach tog $togetherList {
    puts $out "$tog"
}
close $out

set nameList {}
set out [open ${outPrefix}_combined.txt w]
foreach comb $combinedList {
    puts $out "$comb"
    lappend nameList [lindex $comb 0]
}
close $out

set i 0
set out [open ${outPrefix}_unique.txt w]
set out1 [open ${outPrefix}_unique.dat w]
foreach comb $uniList {
    puts $out "$comb"
    puts $out1 "$i $comb"
    incr i
}
close $out

# Sort by reverses.
set reverseList {}
foreach comb $uniList {
    set c [string index [lindex $comb 0] 2]
    # Move to a representation where all sequences have a purine in the middle.
    if {[string equal $c t] || [string equal $c c]} {
	set canon [complement [lindex $comb 0]]
    } else {
	set canon [lindex $comb 0]
    }
    set key [string range $canon 1 3]

    lappend reverseList [concat $canon [lrange $comb 1 2] $key]
}

set reverseList [lsort -index 0 $reverseList]
set reverseList [lsort -index 3 $reverseList]
set diffList {}
foreach {x y} $reverseList {
    set d [expr {[lindex $y 1] - [lindex $x 1]}]
    lappend diffList [list [lindex $x 0] [lindex $x 1] [lindex $x 2] $d]
    lappend diffList [list [lindex $y 0] [lindex $y 1] [lindex $y 2] $d]
}

set out [open ${outPrefix}_reverse.txt w]
foreach comb $diffList {
    puts $out "$comb"
}
close $out

set allList {}
foreach p {+ -} {
    foreach a {a t g c} {
	foreach b {a t g c} {
	    foreach c {a t g c} {
		lappend allList ${p}${a}${b}${c}
	    }
	}
    }
}

set out [open ${outPrefix}_left.txt w]
foreach name $allList {
    if {[lsearch $nameList $name] < 0} {
	puts -nonewline $out "$name "
    }
}
close $out
