# Make a sharp grid based on the chosen atoms.
# Use with: tclsh sharpPhantomGrid.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set radius0 1.7; # switching occurs from radius0 to radius0+sigma
set sigma 0.8

set s0 12
set pi [expr 4.0*atan(1.0)]
set slope [expr tan(10.0*$pi/180.0)]

set gridX 0.6
set gridY 0.6
set gridZ 0.6
set lx 72.0
set ly 72.0
set lz 206.0

#Input:
set inFile pore2.0_coords.txt
#Output:
set outFile pore2.0_phtm${sigma}-${radius0}.dx

source gridForce.tcl
source cellDecomposition.tcl
set cutoff [expr $radius0+$sigma+0.5]
set cutoff2 [expr 2.0*$cutoff]

proc getPotential {r a sigma} {
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

# Load the system and extract the coordinates.
foreach i {0} {set pos [loadPoints $inFile]}
set n [llength $pos]
puts "Read $n points."
set bounds [getBounds $pos]

# Initialize the cells.
puts "Computing cell decomposition..."
set l [vecSub [lindex $bounds 1] [lindex $bounds 0]]
set nx [expr int(ceil([lindex $l 0]/$cutoff))+2]
set ny [expr int(ceil([lindex $l 1]/$cutoff))+2]
set nz [expr int(ceil([lindex $l 2]/$cutoff))+2]
set nCell [list $nx $ny $nz]
set oCell [vecSub [lindex $bounds 0] [list $cutoff $cutoff $cutoff]]
foreach i {0} {set neigh [cellNeighbors $nCell]}
cellDecompose $pos $cutoff $oCell $nCell cell

puts "Cell decomposition dimensions: $nCell"
puts "Cell decomposition origin: $oCell"
puts "Cell decomposition contains [cellCountPoints $nCell cell] points"

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
	    set r [vecAdd [vecTransform $delta [list $ix $iy $iz]] $origin]
	    
	    foreach {x y z} $r {break}
	    set s [expr sqrt($x*$x + $y*$y)]
	    set sPore [expr $s0 + $slope*abs($z)]
	    set pot 0.0
	    
	    if {$s < $sPore - $cutoff2} {
		set pot 0.0
	    } elseif {$s > $sPore + $cutoff2 && abs($z-100.0) > $cutoff2 && abs($z+100) > $cutoff2} {
		set pot 0.0
	    } else {
		# Find our cell.
		set home [cellLookup $r $cutoff $oCell $nCell]
	    
		# Loop through the neighbors of our cell.
		foreach c [lindex $neigh $home] {
		    foreach p $cell($c) {
			set dr [vecLength [vecSub $r $p]]
			set pot [expr $pot + [getPotential $dr $radius0 $sigma]]
		    }
		}
	    }

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



