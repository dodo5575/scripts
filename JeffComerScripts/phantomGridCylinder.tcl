# Make a sharp grid based on the chosen atoms.
# Use with: tclsh sharpPhantomGrid.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set radius0 20.0; # switching occurs from radius0 to radius0+sigma
set sigma 4.0

set radiusPore 10.0
set lengthPore 240.0

set pi [expr 4.0*atan(1.0)]
set slope [expr tan(10.0*$pi/180.0)]

set gridX 2.5
set gridY 2.5
set gridZ 2.5
set lx 140.0
set ly 140.0
set lz 260.0

#Output:
set outFile cylinder${radiusPore}_low.dx

source gridForce.tcl
set s0 [expr $radiusPore + 0.5*$sigma]
set z0 [expr 0.5*$lengthPore - 0.5*$sigma]
set rc $radius0

proc getPotential {r a sigma} {
    if {[expr $r < 0]} {return [expr $a+0.5*$sigma]}
    
    if {[expr $r < $a]} {
	set pot [expr $a - $r + 0.5*$sigma]
    } elseif {[expr $r > $a + $sigma]} {
	set pot 0.0
    } else {
	set dr [expr $r-$a]
	set pot [expr $dr*(0.5*$dr/$sigma-1.0) + 0.5*$sigma]
    }
    
    return $pot
}

# Initialize the grid.
set grid(data) {}
set grid(origin) [vecScale -0.5 [list $lx $ly $lz]]
set grid(nx) [expr int(ceil($lx/$gridX))]
set grid(ny) [expr int(ceil($ly/$gridY))]
set grid(nz) [expr int(ceil($lz/$gridZ))]
set grid(size) [expr $grid(nx)*$grid(ny)*$grid(nz)]
set grid(delta) [list [list $gridX 0 0] [list 0 $gridY 0] [list 0 0 $gridZ]]
set grid(deltaInv) [matInvert $grid(delta)]

puts "Grid dimensions: $grid(nx) $grid(ny) $grid(nz)"
puts "Grid nodes: $grid(size)"
puts "Grid origin: $grid(origin)"

# Write the grid header.
set out [open $outFile w]
writeDxHeader grid $out

# Sample the potential.
puts "Sampling the potential..."
set delta $grid(delta)
set origin $grid(origin)
set complete 0
set j 0
for {set ix 0} {$ix < $grid(nx)} {incr ix} {
    for {set iy 0} {$iy < $grid(ny)} {incr iy} {
	for {set iz 0} {$iz < $grid(nz)} {incr iz} {
	    set p [vecAdd [vecTransform $delta [list $ix $iy $iz]] $origin]
	    foreach {x y z} $p {break}
	    set s [expr sqrt($x*$x + $y*$y)]

	    if {[expr ($s < $s0 - $sigma) || ($z > $z0 + $sigma) || ($z < -$z0 - $sigma)]} {
		set pot 0.0
	    } else {
		if {[expr $s < $s0+$rc && $z > $z0-$rc]} {
		    # Top corner.
		    set ds [expr $s - ($s0+$rc)]
		    set dz [expr $z - ($z0-$rc)]
		    set r [expr sqrt($ds*$ds + $dz*$dz)]
		} elseif {[expr $s < $s0+$rc && $z < -$z0+$rc]} {
		    # Bottom corner.
		    set ds [expr $s - ($s0+$rc)]
		    set dz [expr $z - (-$z0+$rc)]
		    set r [expr sqrt($ds*$ds + $dz*$dz)]
		} elseif {[expr $z < $z0 - $rc && $z > -$z0 + $rc]} {
		    # Inside the pore.
		    set r [expr ($s0+$rc)-$s]
		} else {
		    # Above or below.
		    if {$z > 0} {
			set r [expr $z-($z0-$rc)]
		    } else {
			set r [expr (-$z0+$rc)-$z]
		    }
		}
		
		set pot [getPotential $r $rc $sigma]
	    }
	    
	    # Write the potential.
	    puts -nonewline $out $pot
	    if {[expr $j % 3 == 2]} {
		puts $out ""
	    } else {
		puts -nonewline $out " "
	    }
	    incr j
	}
	
	if {[expr 100*$j/$grid(size) - $complete] >= 5} {
	    set complete [expr 100*$j/$grid(size)]
	    puts "Percent complete: $complete"
	}
    }
}

# Write the result.
close $out
puts "$j potential values were written."
exit



