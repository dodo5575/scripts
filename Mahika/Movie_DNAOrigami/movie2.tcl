## User-defined movie frame update callback procedure, invoked every time
## a movie frame is rendered.

set trajectory_frames [molinfo top get numframes]

proc moviecallback { args } {
    global trajectory_frames
    puts "User-defined movie frame update callback frame: $::MovieMaker::userframe
 / $::MovieMaker::numframes"


    #before running this script, hit '=' while in the graphical display and run setting2.tcl


    animate goto $::MovieMaker::userframe;

    if {$::MovieMaker::userframe <= 75} {
        scale by 1.02
    } else {
        rotate z by 0.2;
        rotate y by 0.2;
    }

}

## Easy-to-use proc to enable the user-defined movie frame callback
proc enablemoviecallback { }  {
    animate goto 0
    trace add variable ::MovieMaker::userframe write moviecallback
}

## Easy-to-use proc to disable the user-defined movie frame callback
proc disablemoviecallback { } {
  trace remove variable ::MovieMaker::userframe write moviecallback
}

