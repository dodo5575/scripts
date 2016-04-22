## User-defined movie frame update callback procedure, invoked every time
## a movie frame is rendered.

set trajectory_frames [molinfo top get numframes]

proc moviecallback { args } {
    global trajectory_frames
    puts "User-defined movie frame update callback frame: $::MovieMaker::userframe
 / $::MovieMaker::numframes"

    if {$::MovieMaker::userframe >= 0} {   
        rotate x by 0.2;
    }

    animate goto $::MovieMaker::userframe;
    

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

