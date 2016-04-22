# A lot of procs that I use a lot.
# Author: Jeff Comer <jeffcomer at gmail>
source $env(HOME)/scripts/vector.tcl

proc trimPath {name} {
    set ind [string last "/" $name]
    return [string range $name [expr {$ind+1}] end]
}

proc extractPath {name} {
    set ind [string last "/" $name]
    return [string range $name 0 [expr {$ind-1}]]
}

proc trimExtension {name} {
    set ind [string last "." $name]
    return [string range $name 0 [expr {$ind-1}]]
}

proc silent {} {
}

# Just read a space delimited data file.
proc readData {fileName} {
    set in [open $fileName r]
    
    set r {}
    while {[gets $in line] >= 0} {
	if {[string match "#*" $line]} { continue }
	if {[string length [string trim $line]] < 1} { continue }

	set tok [concat $line]
	lappend r $tok
    }

    close $in
    return $r
}

proc writeData {fileName data} {
    set out [open $fileName w]
    
    foreach d $data {
	puts $out $d
    }

    close $out
}

# Construct a pdb line from everything.
proc makePdbLineBeta {index segName resId name resName r beta} {
    set template "ATOM     42  HN  ASN X   4     -41.083  17.391  50.684  0.00  0.00      P1    "

    foreach {x y z} $r {break}
    set record "ATOM  "
    set si [string range [format "     %5i " $index] end-5 end]
    if {[string length $name] < 4} {
	set name [string range " $name    " 0 3]
    } else {
	set name [string range $name 0 3]
    }
    set resName [string range " $resName     " 0 4]
    set chain "[string index $segName 0]"
    set resId [string range "    $resId"  end-3 end]
    set temp1 [string range $template  26 29]
    set sx [string range [format "       %8.3f" $x] end-7 end]
    set sy [string range [format "       %8.3f" $y] end-7 end]
    set sz [string range [format "       %8.3f" $z] end-7 end]
    set temp2 [string range $template 54 59]
    set beta [string range [format "       %6.2f" $beta] end-5 end]
    set temp3 [string range $template 66 71]
    set segName [string range "$segName    "  0 3]
    set tempEnd [string range $template 76 end]

    # Construct the pdb line.
    return "${record}${si}${name}${resName}${chain}${resId}${temp1}${sx}${sy}${sz}${temp2}${beta}${temp3}${segName}${tempEnd}"
}

# Dave's vmdargs
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


proc ul { str } {
    # This proc brackets the provided string with the terminal escape sequence for underlining
    return "\033\[4m$str\033\[0m"
}

### STACK PROCS ###
# from Dave
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

