# Author: Jeff Comer <jcomer2@illinois.edu>
set dcdDir /scratch/tbgl/jcomer/basepair/dcd
set numInstances 8
set bp null
set field pos
set pattern "b*${bp}_${field}*.dcd"
set outName dcdList_${bp}_${field}

# Get rid of nonmatching list elements.
proc lnmatch {lst reg} {
    set ret {}
    foreach item $lst {
	if {![regexp $reg $item]} {lappend ret $item}
    }
    return $ret
}

proc writeList {name lis} {
     set out [open $name w]
    foreach item $lis {
	puts $out $item
    }
    close $out
}


# Glob the dcd files, skipping simulations possibly not in a steady state.
set dcdList [glob -directory $dcdDir $pattern]
#set dcdList [lnmatch $dcdList  ".*basepair_.*"]
#set dcdList [lnmatch $dcdList  ".*bp0_.*"]
set dcdList [lsort $dcdList]
puts "There are a total of [llength $dcdList] files."

set n [expr int([llength $dcdList]/$numInstances)]
for {set i 0} {$i < [expr $numInstances-1]} {incr i} {
    set j0 [expr $n*$i]
    set j1 [expr $n*($i+1) - 1]
    
    set lis [lrange $dcdList $j0 $j1]
    writeList ${outName}${i}.txt $lis
}
set lis [lrange $dcdList [expr $n*($numInstances-1)] end]
writeList ${outName}${i}.txt $lis
