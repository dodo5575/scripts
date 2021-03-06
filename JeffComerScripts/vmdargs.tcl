proc ul { str } {
    # This proc brackets the provided string with the terminal escape sequence for underlining
    return "\033\[4m$str\033\[0m"
}

proc vmdargs { script args } {
    upvar argv argv_local
    
    if { [llength $argv_local] != [llength $args] } {
	puts -nonewline "Usage: vmd -dispdev text -e $script -args"
	foreach arg $args {
	    puts -nonewline " [ul $arg]"
	}
	puts ""
	exit -1
    }
    
    foreach arg $args {
	upvar $arg arg_local
	set arg_local [shift argv_local]
    }
}

proc vmdargslist { script args } {
    upvar argv argv_local
    
    # DEBUG
    puts "argv: $argv_local"
    
    set nargs [llength $args]
    set nlistargs 0
    set i 0
    foreach arg $args {
	# Loop through and check for a list arg ... more than one is a no-no
	if { [string match "@*" $arg] } {
	    incr nlistargs
	    set listargind $i
	}
	incr i
    }
    if { $nlistargs > 1 } {
	puts "Cannot use more than one list variable in vmdargs call!"
	exit -1
    }
    
    if { ($nlistargs == 0 && [llength $argv_local] != [llength $args]) || [llength $argv_local] < [llength $args] } {
	puts -nonewline "Usage: vmd -dispdev text -e $script -args"
	foreach arg $args {
	    if { [string match "@*" $arg] } {
		set arg2 [string range $arg 1 end]
		puts -nonewline " [ul $arg2] \[ [ul $arg2] ... \]"
	    } else {
		puts -nonewline " [ul $arg]"
	    }
	}
	puts ""
	exit -1
    }
    
    foreach arg $args {
	if { [string match "@*" $arg] } {
	    # remove leading '@' symbol
	    set arg2 [string range $arg 1 end]
	    upvar $arg2 arg_local
	    while { [llength $argv_local] >= $nargs - $listargind } {
		lappend arg_local [shift argv_local]
	    }
	} else {
	    upvar $arg arg_local
	    set arg_local [shift argv_local]
	}
    }
}


### STACK PROCS ###
#
proc push { stack value } {
    upvar $stack list
    lappend list $value
}

proc pop { stack } {
    upvar $stack list
    set value [lindex $list end]
    set list [lrange $list 0 [expr [llength $list]-2]]
    return $value
}

proc shift { stack } {
    upvar $stack list
    set value [lindex $list 0]
    set list [lrange $list 1 end]
    return $value
}

proc unshift { stack value } {
    upvar $stack list
    set list [concat $value $list]
}

proc stackrotate { stack {value 1} } {
    upvar $stack list
    while { $value > 0 } {
	set el [shift list]
	push list $el
	incr value -1
    }
    while { $value < 0 } {
	set el [pop list]
	unshift list $el
	incr value 1
    }
    return $list
}
