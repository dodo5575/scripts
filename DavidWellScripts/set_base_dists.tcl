proc enablebasedist {} {
    global vmd_frame
    global sels resids nresids
    
#     set sel_dna [atomselect top "DNA"]
#     set resids [lsort -integer -unique [$sel_dna get resid]]
#     set nresids [llength $resids]
#     $sel_dna delete
#     set sels {}
#     foreach resid $resids {
# 	set sel [atomselect top "DNA and resid $resid"]
# 	lappend sels $sel
#     }
    
    trace variable vmd_frame([molinfo top]) w drawcounter
}

proc disablebasedist {} {
    global vmd_frame
    global sels resids nresids
    
    foreach sel $sels {
	$sel delete
    }
    
    trace vdelete vmd_frame([molinfo top]) w drawcounter
}

set sel_dna [atomselect top "DNA"]
set resids [lsort -integer -unique [$sel_dna get resid]]
set nresids [llength $resids]
$sel_dna delete
set sels {}
foreach resid $resids {
    set sel [atomselect top "DNA and resid $resid"]
    lappend sels $sel
}

proc drawcounter { name element op } {
    global vmd_frame
    global sels resids nresids
    
    set xs {}
    foreach sel $sels {
	lappend xs [measure center $sel weight mass]
    }
    
    set diffs {}
    foreach x $xs {
	if { [info exists last_x] } {
	    lappend diffs [veclength [vecsub $x $last_x]]
	}
	set last_x $x
    }
    
    set i 0
    foreach sel $sels {
	if { $i == 0 } {
	    $sel set user [lindex $diffs 0]
	} elseif { $i == [expr $nresids-1] } {
	    $sel set user [lindex $diffs [expr $nresids-2]]
	} else {
	    $sel set user [expr ([lindex $diffs [expr $i-1]] + [lindex $diffs $i]) / 2]
	}
	incr i
    }
}
