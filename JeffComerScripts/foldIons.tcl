# Get the system size using measure minmax.
# to use: vmd -dispdev text -e getSize.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set x0 -50
set x1 50
set y0 -50
set y1 50
set z0 55
set z1 180
set selText "ions and not within 3.0 of water"
set goodText "not within 2.0 of not (water or ions)"
set margin 0.5
# Input:
set psf eco_start_sys.psf
set pdb eco_start_sys.pdb
# Output:
set outPdb eco_start_sys_fold.pdb

mol load psf $psf pdb $pdb
set sel [atomselect top $selText]
foreach zero {0} {set indexList [$sel get index]}
set nSel [$sel num]
$sel delete

puts "Folding $nSel ions."
set n 0
set reportPeriod 20

foreach i $indexList {
    set s [atomselect top "index $i"]
    set goodNum 0
    
    if {$n % $reportPeriod == 0} {
	puts [format "%.3f percent complete" [expr 100.0*$n/$nSel]]
    }

    while {!$goodNum} {
	set x [expr $x0 + rand()*($x1-$x0)]
	set y [expr $y0 + rand()*($y1-$y0)]
	set z [expr $z0 + rand()*($z1-$z0)]
	$s set {x y z} [list [list $x $y $z]]
	
	set sv [atomselect top "index $i and ($goodText) and not within $margin of (not index $i)"]
	set goodNum [$sv num]
	$sv delete
    }

    $s delete
    incr n
}

set all [atomselect top all]
$all writepdb $outPdb
$all delete

mol delete top
exit



