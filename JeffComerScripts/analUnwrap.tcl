# Author: Jeff Comer <jcomer2@illinois.edu>
# Unwrap given some files.
if {$argc < 5} {
    puts "$argv0 name structPrefix outputDir stride dcdFile0 [dcdFile1...]"
    exit
}
package require pbctools

set name [lindex $argv 0]
set struct [lindex $argv 1]
set outDir [lindex $argv 2]
set stride [lindex $argv 3]
set dcdList [lrange $argv 4 end]

proc compute {name struct outDir stride dcdList} {
    set psf $struct.psf

    mol load psf $psf
    foreach dcd $dcdList {
	mol addfile $dcd step $stride waitfor all
    }

    set nFrames [molinfo top get numframes]
    puts [format "Reading %i frames." $nFrames]
    set last [expr {$nFrames-1}]

    pbc unwrap -all -sel "segname \"L.*\""
    puts "Unwrapped."

    set all [atomselect top all]
    animate write dcd unwrap${stride}_${name}.dcd beg 0 end $last waitfor all sel $all top
    $all delete
    mol delete top
}

compute $name $struct $outDir $stride $dcdList
exit
