# Author: Jeff Comer <jcomer2@illinois.edu>
# Distribute ions randomly inside membrane.

mol delete all
# Parameters:
set poreLength 100.0
set poreDiameter 50.0
set poreAngle 10.0
set basisFile ../0_expel_grid/grid_basis.txt
set psf sio.psf
set pdb sio.pdb
set outPdb sio_ready.pdb

set pi [expr {4.0*atan(1.0)}]
set halfLen [expr {0.5*$poreLength}]
set s0 [expr {0.5*$poreDiameter}]
set slope [expr {tan($poreAngle*$pi/180.0)}]

proc inPoreWalls {pos} {
    global halfLen s0 slope
    
    foreach {x y z} $pos { break }
    if {abs($z) > $halfLen} { return 0 }
    
    set s [expr {sqrt($x*$x + $y*$y)}]
    if {$s < $s0 + $slope*abs($z)} { return 0 }

    return 1
}

# Load the basis vectors.
set cellVecList {}
set in [open $basisFile r]
foreach line [split [read $in] "\n"] {
    lappend cellVecList [concat $line]
    puts "cellBasisVector $line"
}
close $in
set cellA [lindex $cellVecList 0]
set cellB [lindex $cellVecList 1]
set cellC [lindex $cellVecList 2]

# Load the molecule.
mol load psf $psf pdb $pdb
set all [atomselect top all]
set num [$all num]

# Reposition the atoms.
puts "Repositioning atoms."
set posList {}
for {set i 0} {$i < $num} {incr i} {
    set bad 1

    # Find positions within the membrane.
    while {$bad} {
	set a [vecscale [expr {rand()-0.5}] $cellA]
	set b [vecscale [expr {rand()-0.5}] $cellB]
	set c [vecscale [expr {rand()-0.5}] $cellC]
	set r [vecadd $a $b $c]

	if {[inPoreWalls $r]} {
	    set bad 0
	}
    }
    lappend posList $r
}

# Set the positions and write the pdb.
$all set {x y z} $posList
$all writepdb $outPdb

exit
