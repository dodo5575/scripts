# Author: Jeff Comer <jcomer2@illinois.edu>
set stride 10
set nFrames [molinfo top get numframes]
set prefix frames/strept_movie
set count 0

proc renderFrameTachyon {prefix} {
    render Tachyon $prefix.dat ""
}

proc renderFrameSnapshot {prefix} {
    render snapshot $prefix.tga ""
}

proc renderFrame {} {
    global count prefix

    display update
    renderFrameSnapshot ${prefix}${count}
    #renderFrameTachyon ${prefix}${count}
    puts "frame $count"
    
    incr count
    return $count
}

for {set f 0} {$f < $nFrames} {incr f $stride} {
    molinfo top set frame $f
    renderFrame
}
