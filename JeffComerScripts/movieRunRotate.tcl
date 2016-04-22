##
## A very simple example of a user-defined movie script
##
## This example rotates a structure, and fades it in as the movie
## progresses.
##
## To try it out, run the "setupmovie" proc (see the end of this script)
## which loads a structure, displays it with VDW rep, creates and assigns
## a new material, and enables the user-defined movie callback.
## Once that's ready, you can open the movie maker plugin, set the 
## movie type to user-defined, and click go and you should see the script
## put through its' paces.
set f 0.0
set rate 0.05
set count 0

# fade in a material as the movie progresses
proc dorotate {} {
   	rotate x by 0.01
	rotate y by 0.2
}

## update the display as the movie runs
proc moviecallback { args } {
global f
global rate
global count

  puts "User-defined movie frame update callback frame: $::MovieMaker::userframe
 / $::MovieMaker::numframes"
  
  incr count
  if { $count < 100 } {
  	dorotate
	return
	}
  
  set nFrames [molinfo top get numframes]
  set f [expr $f+$rate]
  set curr [expr int($f)]
  
  if {$curr < $nFrames} {
  	molinfo top set frame $f
	} else {
	dorotate
	}
  
}

## Easy-to-use proc to enable the user-defined movie frame callback
proc enablemoviecallback { }  {
  trace add variable ::MovieMaker::userframe write moviecallback
}

## Easy-to-use proc to disable the user-defined movie frame callback
proc disablemoviecallback { } {
  trace remove variable ::MovieMaker::userframe write moviecallback
}

## setup for making the movie
proc setupmovie { } {
  global f
  global count
  
  set f 0.0
  set count 0
  enablemoviecallback
  after idle { puts "Don't forget to set movie type to user-defined in the Movie Settings" }
}  

## disable movie callback, any other cleanup....
proc finishmovie {} {
  disablemoviecallback
}










