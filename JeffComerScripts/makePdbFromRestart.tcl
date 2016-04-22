# Author: Jeff Comer <jcomer2@illinois.edu>
set prefix flow
set vel {10 100 500}
set surf {raw middling anneal}
set conc {10mM 100mM}
set restartList {}
foreach v $vel {
    foreach s $surf {
	foreach c $conc {
	    set f ${prefix}${v}_${s}_${c}.1.restart

	    if {[file exists $f]} { lappend restartList $f }
	}
    }
}

# Input
set outDir poise
set stride 1
set startFrame 0

# Negative frames count from the end.
#set frameList {-1}

proc trimExtension {name} {
    set ind [string last "." $name]
    return [string range $name 0 [expr $ind-1]]
}

# Construct a pdb line from everything.
proc makePdbLineFull {index segName resId name resName r} {
    set template "ATOM     42  HN  ASN X   4     -41.083  17.391  50.684  0.00  0.00      P1    "

    foreach {x y z} $r {break}
    set record "ATOM  "
    set si [string range [format "     %5i " $index] end-5 end]
    if {[string length $name] < 4} {
	set name [string range " $name    " 0 3]
    } else {
	set name [string range $name 0 3]
    }
    set resName [string range " $resName    " 0 3]
    set temp0 " [string index $segName 0]"
    set resId [string range "    $resId"  end-3 end]
    set temp1 [string range $template  26 29]
    set sx [string range [format "       %8.3f" $x] end-7 end]
    set sy [string range [format "       %8.3f" $y] end-7 end]
    set sz [string range [format "       %8.3f" $z] end-7 end]
    set temp2 [string range $template 54 71]
    set segName [string range "$segName    "  0 3]
    set tempEnd [string range $template 76 end]

    # Construct the pdb line.
    return "${record}${si}${name}${resName}${temp0}${resId}${temp1}${sx}${sy}${sz}${temp2}${segName}${tempEnd}"
}

proc readRestart {fileName} {
    set inp [open $fileName r]

    set ret {}
    while {[gets $inp line] >= 0} {
	if {[string match "#*" $line]} {continue}
	if {[string length $line] < 2} {continue}
	
	set tok [concat $line]
	if {[llength $tok] != 4} {
	    puts "readRestart: Invalid format"
	}
	lappend ret $tok
    }

    close $inp
    return $ret
}

proc makePdb {atomList outName} {
    set out [open $outName w]
    set i 0

    foreach atom $atomList {
	foreach {type x y z} $atom { break }
	set r [list $x $y $z]
	set atomId [getAtomId $i $type]
	foreach {segName resId name} $atomId { break }

	set line [makePdbLineFull $i $segName $resId $name PAR $r]  
	puts $out $line
	incr i
    }
    close $out

    return
}

proc getAtomId {index type} {
    set nameList {C O N P S}
    
    set nameInd [expr {$type % [llength $nameList]}]
    set name [lindex $nameList $nameInd]
    set resId [expr {$index % 1000}]
    set segNum [expr {$index/1000}]
    
    if {$segNum >= 1000} {
	puts "ERROR: Too many segment names!"
	exit
    }
    set segName "P${segNum}"
    
    return [list $segName $resId $name]
}


foreach f $restartList {
    set atomList [readRestart $f]
#    puts "length: [llength $atomList]"
#    puts "length: [llength [lindex $atomList 0]]"
    makePdb $atomList $outDir/$f.pdb
}

exit
