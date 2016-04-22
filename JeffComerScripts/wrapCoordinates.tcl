# Wrap coordinates to a basis.
# Use with: vmd -dispdev text -e wrapCoordinates.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set a {60.76 0.0 0.0}
set b {30.38 52.6197035339 0.0}
set c {0.0 0.0 200.0}
#Input:
set psf pore2.0_long.psf
set coords pore2.0_long.pdb
# Output:
set outFile pore2.0_long_wrap.pdb

source vector.tcl 

mol load psf $psf
mol addfile $coords
set sel [atomselect top all]
foreach zero {0} {set pos [$sel get {x y z}]}
set n [$sel num]
puts "\nWrapping $n atoms..."

# Compute the origin.
set r0 [vecAdd $a $b]
set r0 [vecAdd $r0 $c]
set origin [vecScale -0.5 $r0]

# Compute the basis.
set delta [matTranspose [list $a $b $c]]
set deltaInv [matInvert $delta]

set pw {}
foreach r $pos {
    # Transform into lattice space.
    set l [vecTransform $deltaInv [vecSub $r $origin]]
    foreach {lx ly lz} $l {break}
    
    # Wrap.
    set lx [expr $lx - int(floor($lx))]
    set ly [expr $ly - int(floor($ly))]
    set lz [expr $lz - int(floor($lz))]

    # Transform back into world space.
    set r1 [vecAdd [vecTransform $delta [list $lx $ly $lz]] $origin]
    
    lappend pw $r1
}

$sel set {x y z} $pw

$sel writepdb $outFile
$sel delete
mol delete top
exit



