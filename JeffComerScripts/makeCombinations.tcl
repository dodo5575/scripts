# Author: Jeff Comer <jcomer2@illinois.edu>
set ion pot
set baseGrid base_zero.dx
set templateName "ade"
# Input:
set inTrans struct_single_aaa_${ion}.txt
# Output:
set outPrefix combo/single

set baseList {A T G C M}
set baseName(A) ade
set baseName(T) thy
set baseName(G) gua
set baseName(C) cyt
set baseName(M) mcyt

source $env(HOME)/scripts/useful.tcl

set seqList {}
foreach b0 $baseList {
    foreach b1 $baseList {
	foreach b2 $baseList {
	    lappend seqList "${b0}${b1}${b2}"
	}
    }
}

set data [readData $inTrans]

foreach seq $seqList {
    # Open the output transform file and write the first line.
    set out [open ${outPrefix}_${seq}_${ion}.trans w]
    puts $out "$baseGrid 1.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0"

    puts $seq
    for {set i 1} {$i <= 3} {incr i} {
	# Get the base name for this sequence.
	set thisName $baseName([string index $seq [expr {$i-1}]])

	# Extract the data from this line of the file.
	set trans [lrange [lindex $data $i] 1 end]
	set inGrid [lindex [lindex $data $i] 0]
	
	# Substitute the base name.
	set outGrid [regsub $templateName $inGrid $thisName]
	
	# Write the new line.
	puts $out "$outGrid $trans"
    }

    close $out
}
