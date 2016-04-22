# Author: Jeff Comer <jcomer2@illinois.edu>
set psf nw_trap2.0_0.1M.psf
set pdb0 last.pdb
set pdb1 shiftAB.pdb
set numFrames 20
set selText "segname ADNA BDNA"
set outPrefix step_frames/step

mol load psf $psf pdb $pdb0
set sel [atomselect top $selText]
set all [atomselect top all]
foreach quiet {0} { set pos0 [$sel get {x y z}] }
mol addfile $pdb1
foreach quiet {0} { set pos1 [$sel get {x y z}] }

for {set f 0} {$f < $numFrames} {incr f} {
    set x1 [expr {double($f)/$numFrames}]
    set x0 [expr {1.0-$x1}]

    set pos {}
    foreach r0 $pos0 r1 $pos1 {
	set r [vecadd [vecscale $x0 $r0] [vecscale $x1 $r1]]
	#set r [list [lindex $r0 0] [lindex $r0 1] [expr {[lindex $r0 2]+$f}]]
	#set r [vecscale $x0 $r0]
	lappend pos $r
    }
    $sel set {x y z} $pos
    $all writepdb ${outPrefix}${f}.pdb
}

$sel delete
$all delete
mol delete top
exit

