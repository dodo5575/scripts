#!/usr/bin/tclsh

set symList [list A T C G]

set comb {}
foreach sym0 $symList {
    foreach sym1 $symList {
	foreach sym2 $symList {
	    lappend comb ${sym0}${sym1}${sym2}
	}
    }
}

proc transform {s} {
    set ret [string map {T A A T C G G C} $s]
    set ret [string reverse $s]
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

# Complement a sequence.
proc complement {s} {
    set n [string length $s]
    set ret ""

    for {set i 0} {$i < $n} {incr i} {
	set c [string index $s $i]
	
	if {[string equal $c A]} {
	    set c1 T
	} elseif {[string equal $c T]} {
	    set c1 A
	} elseif {[string equal $c G]} {
	    set c1 C
	} elseif {[string equal $c C]} {
	    set c1 G
	} else {
	    puts "Warning: Unrecognized nucleotide $c!"
	    set c1 $c
	}
	
	set ret ${c1}${ret}
    }
    return $ret
}

proc compPairs {combList} {
    set n [llength $combList]
    set remList {}

    for {set i 0} {$i < [expr $n-1]} {incr i} {
	for {set j [expr $i + 1]} {$j < $n} {incr j} {
	    set si [lindex $combList $i]
	    set sj [lindex $combList $j]
	    set tsj [complement $sj]
	    set rsj [reverse $sj]
	    if {[string equal $si $sj] || [string equal $si $tsj] || [string equal $si $rsj]} {
		lappend remList $si
	    } 
	}
    }
    
    set retList {}
    foreach s $combList {
	if {[lsearch $remList $s] < 0} {
	    lappend retList $s
	}
    }
    return $retList
}


set len 0
while {[llength $comb] != $len} {
    set len [llength $comb]
    set comb [compPairs $comb]
    puts $len
}

puts $comb


