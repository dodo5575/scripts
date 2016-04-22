#
# Auth: RCarr 11/19/2008
#
# Draws WHAM windows.
#

source /home/rccarr2/Scripts/Common.tcl



proc read_windows { inFile {temperature 295} } {
    global Windows Coverage kB_JperK JtoKCalperMol Box

    # Now, get the forces
    set FH [open $inFile r]
    while {[gets $FH line] >= 0} {
	if {[string match "\#*" $line]} {
	    if { [llength $line] == 7 } {
		set Box [lrange $line 1 6]
		puts "Using Box $Box"
	    }
	    continue 
	}

	if {[llength $line] < 7} { continue }
			
	set Windows([lindex $line 0]) [lrange $line 1 3]
#	set Springs([lindex $line 0]) [lrange $line 4 6]

	# only supports spherical windows right now.
	set Coverage([lindex $line 0]) [expr { sqrt($kB_JperK * $JtoKCalperMol*$temperature/[lindex $line 4]) }]
	
    }
    close $FH

    puts "Read in [array size Windows] windows"

}

proc draw_windows { {radius 0.2} } {
    global Windows MolWindows
    lappend MolWindows [mol new]
    foreach name [array names Windows] {
	draw sphere $Windows($name) resolution 27 radius $radius
    }
}

proc draw_coverage { {scale 2.0} } {
    global Windows Coverage MolWindows DrawRef
    lappend MolWindows [mol new]
    draw color orange
    foreach name [array names Windows] {
	set i [draw sphere $Windows($name) resolution 27 radius [expr {$scale * $Coverage($name)}]]
	set DrawRef($name) $i
    }
}

proc draw_box { {thickness 3} } {
    global Box MolWindows

    lappend MolWindows [mol new]

    set origin [lrange $Box 3 5]

    draw line $origin [vecadd $origin [list 0.0 0.0 [lindex $Box 2]] ] width $thickness
    draw line $origin [vecadd $origin [list 0.0 [lindex $Box 1] 0.0] ] width $thickness
    draw line $origin [vecadd $origin [list [lindex $Box 0] 0.0 0.0] ] width $thickness
    
    set basevec [vecadd $origin [list 0.0 0.0 [lindex $Box 2] ] ]
    draw line $basevec [vecadd $basevec [list 0.0 [lindex $Box 1] 0.0] ] width $thickness
    draw line $basevec [vecadd $basevec [list [lindex $Box 0] 0.0 0.0] ] width $thickness

    set basevec [vecadd $origin [list [lindex $Box 0] 0.0 0.0] ]
    draw line $basevec [vecadd $basevec [list 0.0 [lindex $Box 1] 0.0] ] width $thickness
    draw line $basevec [vecadd $basevec [list 0.0 0.0 [lindex $Box 2]] ] width $thickness

    set basevec [vecadd $origin [list [lindex $Box 0] 0.0 [lindex $Box 2] ] ]
    draw line $basevec [vecadd $basevec [list 0.0 [lindex $Box 1] 0.0] ] width $thickness

    set basevec [vecadd $origin [list [lindex $Box 0] [lindex $Box 1] 0.0 ] ]
    draw line $basevec [vecadd $basevec [list 0.0  0.0 [lindex $Box 2]] ] width $thickness

    set basevec [vecadd $origin [list  0.0 [lindex $Box 1] [lindex $Box 2] ] ]
    draw line $basevec [vecadd $basevec [list [lindex $Box 0] 0.0  0.0] ] width $thickness    

    set basevec [vecadd $origin [list 0.0 [lindex $Box 1] 0.0] ]
    draw line $basevec [vecadd $basevec [list [lindex $Box 0] 0.0 0.0] ] width $thickness
    draw line $basevec [vecadd $basevec [list 0.0 0.0 [lindex $Box 2]] ] width $thickness
    
    set basevec [vecadd $origin [list 0.0 [lindex $Box 1] 0.0] ]
    draw line $basevec [vecadd $basevec [list [lindex $Box 0] 0.0 0.0] ] width $thickness
    draw line $basevec [vecadd $basevec [list 0.0 0.0 [lindex $Box 2]] ] width $thickness
    

}

proc clear_window_drawings { } {
    global MolWindows
    foreach mol $MolWindows {
	mol delete $mol
    }
}

proc erase_point { ref } {
    global DrawRef
    draw delete $DrawRef($ref)
}

proc reset_windows { } {
    global Windows Coverage DrawRef Box MolWindows

    array unset Windows
    array unset Coverage 
    array unset DrawRef
    
    array set Windows [list ]
    array set Coverage [list ]
    array set DrawRef [list ]
    set Box [list ]
    set MolWindows [list ]

    clear_window_drawings
    
}

reset_windows
