source $env(SCRIPTS)/Procs.tcl

namespace eval ::Dave::dirstack {
    variable dirstack [pwd]
    
    proc pushd { args } { 
	variable dirstack
	
	set value 0
	set dir ""
	if { [llength $args] == 1 } {
	    set arg [lindex $args 0]
	    if { [string is integer -strict $arg] } {
		set value $arg
	    } else {
		set dir $arg
	    }
	}
	if { $dir eq "" && [llength $dirstack] < 2 } {
	    puts "pushd: no other directory"
	    return
	}
	if { [string range $value 0 0] eq "-" } {
	    # need this since $value could be "-0"
	    incr value -1
	}
	
	set value [expr $value % [llength $dirstack]]
	if { $value == 0 } {
	    if { $dir eq "" } {
		set dirstack [lreplace $dirstack 0 1 [lindex $dirstack 1] [lindex $dirstack 0]]	;# swap first and second elements
	    } else {
		unshift dirstack $dir
	    }
	} else {
	    stackrotate dirstack $value
	}
	cd [lindex $dirstack 0]
	
	return $dirstack
    }
    
    proc popd { {value 0} } {
	variable dirstack
	
	if { [llength $dirstack] < 2 } {
	    puts "popd: directory stack empty"
	    return
	}
	
	if { [string range $value 0 0] eq "-" } {
	    # need this since $value could be "-0"
	    incr value -1
	}
	set value [expr $value % [llength $dirstack]]
	if { $value == 0 } {
	    cd [lindex $dirstack 1]
	}
	set dirstack [lreplace $dirstack $value $value]	;# deletes $value-th item
	return $dirstack
    }
    
    proc dirs { } { 
	variable dirstack
	return $dirstack
    }
    
    proc cd_hook { args } {
	variable dirstack
	shift dirstack
	unshift dirstack [pwd]
	#puts $args
    }
    
    namespace export pushd
    namespace export popd
    namespace export dirs
}

namespace import ::Dave::dirstack::*

trace add execution cd leave ::Dave::dirstack::cd_hook
