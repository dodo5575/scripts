# Author: Jeff Comer <jcomer2@illinois.edu>
set stride 1
set nFrames [molinfo top get numframes]
set prefix frames/step_movie
#set colorList {4 20 0 11}
set colorList {30 0 30 27}

set runFrames 145
set stepFrames 21

source representDnaCartoonMod.tcl
set count 0
set nSteps [llength $colorList]

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

set count 0
representDna ADNA BDNA $colorList -1
for {set s 0} {$s < $nSteps} {incr s} {
    
    for {set f 0} {$f < $runFrames} {incr f} {
	molinfo top set frame $f
	renderFrame
    }
    
    for {set f $runFrames} {$f < $nFrames} {incr f} {
	molinfo top set frame $f
	renderFrame
    }
    representDna ADNA BDNA $colorList $s
}
