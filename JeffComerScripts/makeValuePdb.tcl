#!/usr/bin/tclsh
# jcomer2@illinois.edu
if {$argc != 2} {
    puts "$argv0 inFile outFile"
    exit
}
set inFile [lindex $argv 0]
set outFile [lindex $argv 1]

set segName MOL
set name C
set resName MOL

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
    set resName [string range " $resName    " 0 3]
    set temp0 " [string index $segName 0]"
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
    return "${record}${si}${name}${resName}${temp0}${resId}${temp1}${sx}${sy}${sz}${temp2}${beta}${temp3}${segName}${tempEnd}"
}

# Just read a space delimited data file.
proc readData {fileName} {
    set in [open $fileName r]
    
    set r {}
    while {[gets $in line] >= 0} {
	if {[string match "#*" $line]} {continue}
	if {[string length $line] < 2} {continue}

	set tok [concat $line]
	lappend r $tok
    }

    close $in
    return $r
}

set data [readData $inFile]

set out [open $outFile w]
set i 1
foreach field $data {
    set r [lrange $field 0 2]
    set v [lindex $field 3]
    
    set lin [makePdbLineBeta $i $segName $i $name $resName $r $v]
    puts $out $lin
    incr i
}
close $out