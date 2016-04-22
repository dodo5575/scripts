# Author: Jeff Comer <jcomer2@illinois.edu>
set channelLen 200
set channelWidth 100
# Input
set baseGrid small_sys.dx
set barrierGrid cube_barrier_blur.dx
set surfGrid cube_raw_hard.dx
set cornerGrid cube_raw_corner.dx
# Output
set outFile first_sys.def

source $env(HOME)/scripts/gridForce.tcl
source $env(HOME)/scripts/vector.tcl

readDx grid $baseGrid
set origin $grid(origin)
set x [lindex $origin 0]
set dx [lindex $grid(delta) 0 0]
set nx1 [expr {$grid(nx)-1}]
set mx [expr {floor($grid(nx)/2)}]

set my [expr {floor($grid(ny)/2)}]
set dy [lindex $grid(delta) 1 1]

set chanWidth [expr {int(ceil(0.5*$channelWidth/$dy))}]
set chanLen [expr {int(ceil(0.5*$channelLen/$dx))}]

puts "channelWidth/2: $chanWidth nodes"
puts "channelLen/2: $chanLen nodes"
set out [open $outFile w]

# The base grid.
puts $out "-1 $baseGrid 1"

for {set i 0} {$i < $grid(size)} {incr i} {
    foreach {ix iy iz} [latticePosition grid $i] { break }

    # Barriers at ends of system
    if {$ix == 0} {
	# First barrier
	puts $out "$i $barrierGrid +x"
	continue
    } elseif {$ix == $nx1} {
	# Second barrier
	puts $out "$i $barrierGrid -x"
	continue
    }

    # Do the corners.
    if {$ix == $mx - $chanLen && $iy == $my - $chanWidth} {
	puts $out "$i $cornerGrid -x"
	continue
    } elseif {$ix == $mx - $chanLen && $iy == $my + $chanWidth} {
	puts $out "$i $cornerGrid -y"
	continue
    } elseif {$ix == $mx + $chanLen && $iy == $my - $chanWidth} {
	puts $out "$i $cornerGrid +y"
	continue
    } elseif {$ix == $mx + $chanLen && $iy == $my + $chanWidth} {
	puts $out "$i $cornerGrid +x"
	continue
    }
    
    # The membrane surfaces
    if {$iy < $my - $chanWidth || $iy > $my + $chanWidth} {
	if {$ix == $mx - $chanLen} {
	    puts $out "$i $surfGrid -x"
	} elseif {$ix == $mx + $chanLen} {
	    puts $out "$i $surfGrid +x"
	}
    } elseif {$ix > $mx - $chanLen && $ix < $mx + $chanLen} {
	# Inside the channel.
	if {$iy == $my - $chanWidth} {
	    puts $out "$i $surfGrid +y"
	} elseif {$iy == $my + $chanWidth} {
	    puts $out "$i $surfGrid -y"
	}
    }
}
close $out

exit
