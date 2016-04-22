# Author: Jeff Comer <jcomer2@illinois.edu>

source computeBrownTownTraj.tcl

set inDir ..
set inPrefix best

set outDir flux
set analysis computeFlux
set outPrefix flux

#set outDir number
#set analysis computeNumber
#set outPrefix num

proc trimExtension {name} {
    set ind [string last "." $name]
    return [string range $name 0 [expr $ind-1]]
}
proc trimPath {name} {
    set ind [string last "/" $name]
    return [string range $name [expr $ind+1] end]
}

set trajList [glob $inDir/${inPrefix}*traj]
set sim {}
foreach f $trajList {
    set name [trimPath [trimExtension $f]]
    lappend sim [list $name $f]
}

foreach s $sim {
    foreach {name traj} $s { break }

    set outFile $outDir/${outPrefix}_${name}.dat
    #set fileList [glob $inDir/$s.traj]
    #set fileList [concat $fileList [lsort [glob $inDir/$s.*.traj]]]
    set fileList $traj
    if {[llength $fileList] == 0} { continue }
    puts $fileList

    computeTrajList $fileList $analysis $outFile
}
