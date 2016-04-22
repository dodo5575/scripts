#!/usr/bin/tclsh
# Author: Jeff Comer <jcomer2@illinois.edu>

if {[llength $argv] != 4} {
    puts "Usage: $argv0 expandCount psf inDcd outDcd"
    exit
}
set expandCount [lindex $argv 0]
set psf [lindex $argv 1]
set inDcd [lindex $argv 2]
set outDcd [lindex $argv 3]

proc expandDcd {expandCount psf inDcd outDcd} {
    mol load psf $psf
    mol addfile $inDcd waitfor all
    
    set nFrames [molinfo top get numframes]
    puts "Loaded $nFrames frames."

    for {set f 0} {$f < $nFrames} {incr f} {
	molinfo top set frame $f
	for {set i 0} {$i < $expandCount} {incr i} {
	    animate dup top
	}
    }
    puts "Currently [molinfo top get numframes] frames."

    animate delete beg 0 end [expr {$nFrames-1}] top
    puts "After duplication, there are [molinfo top get numframes] frames."

    animate write dcd $outDcd waitfor all top
    puts "Wrote `$outDcd'"

    return
}


expandDcd $expandCount $psf $inDcd $outDcd
exit
