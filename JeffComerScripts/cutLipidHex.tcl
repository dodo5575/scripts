# This script will remove water from psf and pdf outside of a
# hexagonal prism along the z-axis with a vertex along the x-axis.
# use with: vmd -dispdev text -e cutWaterHex.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>
package require psfgen 1.3

set name silica_pore40+dopc
# Input: 
set psf ${name}.psf
set pdb ${name}.pdb 
# Output:
set psfFinal ${name}_hex.psf
set pdbFinal ${name}_hex.pdb

# Parameters:
set maxHeight 1000; # maximum extent in both z directions
set padding 2.7; # 
# from the pore hexagon. The distance is in angstrom
set waterText "segname \"L.*\"";# This is the stuff that is removed.
# This selection forms the basis for the hexagon.
set selText "resname SIO2"
# swapXY = 0: Hexagon vertex along x axis
# swapXY = 1: Hexagon face center along x axis
set swapXY 1

proc quiet {} {}

# Load the molecule.
mol load psf $psf pdb $pdb

# Find the system dimensions.
set sel [atomselect top $selText]
set minmax [measure minmax $sel]
$sel delete

set z0 [expr {$maxHeight + $padding}]
set size [vecsub [lindex $minmax 1] [lindex $minmax 0]]
set sqrt3 [expr sqrt(3.0)]
foreach {sizeX sizeY sizeZ} $size {break}
if {!$swapXY} {
    set x x
    set y y
    foreach {sizeX sizeY sizeZ} $size {break}
} else {
    set x y
    set y x
    foreach {sizeY sizeX sizeZ} $size {break}
}

# This is the hexagon's radius.
set r [expr {0.5*$sizeX + $padding}]
#set r [expr {0.5*$sizeX}]

# Find water outside of the hexagon.
# Check the middle rectangle.
set check "($waterText) and (((abs(${x}) < 0.5*$r and abs(${y}) > 0.5*$sqrt3*$r) or"

# Check the lines forming the nonhorizontal sides.
set check [concat $check "(${y} > $sqrt3*(${x}+$r) or ${y} < $sqrt3*(${x}-$r) or"]
set check [concat $check "${y} > $sqrt3*($r-${x}) or ${y} < $sqrt3*(-${x}-$r)))"]

# Check the height.
set check "$check or abs(z) > $z0)"
set w [atomselect top $check]
set violators [lsort -unique [$w get {segname}]]; quiet
$w delete

# Remove the offending water molecules.
puts "Deleting [llength $violators] offending lipids..."
resetpsf
readpsf $psf
coordpdb $pdb
foreach waterMol $violators {
    delatom [lindex $waterMol 0]
}

writepsf $psfFinal
writepdb $pdbFinal
mol delete top
exit
