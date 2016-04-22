# Author: Jeff Comer <jcomer2@illinois.edu>
set stride 1
set run $dropFrame
set nFrames [molinfo top get numframes]
set prefix frames/trap_tachyon

proc renderFrameTachyon {prefix} {
    render Tachyon $prefix.dat ""
}

proc renderFrameSnapshot {prefix} {
    render snapshot $prefix.tga ""
}

color Display {Background} white
for {set f 0} {$f < $run} {incr f $stride} {
    molinfo top set frame $f

    display update
    renderFrameTachyon ${prefix}${f}
    #renderFrameSnapshot ${prefix}${f}

    puts "frame $f"    
}

color Display {Background} pink
for {set f $run} {$f < $nFrames} {incr f $stride} {
    molinfo top set frame $f

    display update
    renderFrameTachyon ${prefix}${f}
    #renderFrameSnapshot ${prefix}${f}

    puts "frame $f"    
}
