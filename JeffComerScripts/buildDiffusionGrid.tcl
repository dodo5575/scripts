# Author: Jeff Comer <jcomer2@illinois.edu>

# Parameters:
set dl 2.0
set lx 40
set ly 40
set lz 72
set diffusion0 {227 241}
set insideFactor 0.80
set nameList {pot chl}
set rad 10
set z0 18
# Output:
set outName diffusion80_

set diffustion1 {}
foreach d $diffusion0 {
    lappend diffusion1 [expr $insideFactor*$d]
}

source vector.tcl
source gridForce.tcl

foreach name $nameList diff0 $diffusion0 diff1 $diffusion1 {
    newGridFit grid $lx $ly $lz $dl

    set newData {}
    set j 0
    foreach val $grid(data) {
	set r [indexToWorld grid $j]

	foreach {x y z} $r { break }
	set s2 [expr $x*$x + $y*$y]
	if {$s2 < $rad*$rad && abs($z) < $z0} {
	    lappend newData $diff1
	} else {
	    lappend newData $diff0
	}
	incr j
    }

    set grid(data) $newData
    writeDx grid ${outName}${name}.dx
}
