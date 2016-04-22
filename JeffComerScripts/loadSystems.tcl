# Author: Jeff Comer <jcomer2@illinois.edu>

# Input:
#set indexFile index_all.txt
set indexFile master_index.txt
set structSuffix _ions
set outputDir output
set repScript representDnaSpheres1.tcl
set poreScript representPore1.tcl

proc quiet {} {}

# Just read a space delimited data file.
proc readData {fileName} {
    set in [open $fileName r]
    
    set r {}
    while {[gets $in line] >= 0} {
	if {[string match "#*" $line]} {continue}
	if {[string length $line] < 2} {continue}

	set tok [concat $line]
	set t [lindex $tok 0]
	lappend r $tok
    }

    close $in
    return $r
}

set data [readData $indexFile]; quiet
set data [lsort -index 1 $data]; quiet

set lastStruct ""
foreach d $data {
    foreach {name struct init voltage} $d { break }

    if {![string equal $lastStruct $struct]} {
	mol load psf ${struct}${structSuffix}.psf pdb ${struct}${structSuffix}.pdb
	source $poreScript
	molinfo top set active off
	molinfo top set displayed off

	mol load psf ${struct}${structSuffix}.psf
	source $repScript
	molinfo top set active on
	molinfo top set displayed off

	set lastStruct $struct
    }

    mol addfile $outputDir/$name.restart.coor
}
