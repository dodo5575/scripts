# Wrap about orthogonal periodic boundaries.
# vmd -dispdev text -e wrapAtoms.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

# parameters:
set lx 79.648
set ly 79.648
set lz 150.0
set selText all
# input:
set psfIn silica_shell_big.psf
set pdbIn silica_shell_big.pdb
set coords min1.restart.2000.coor
# output:
set pdbOut min1_wrap.pdb

mol load psf $psfIn pdb $pdbIn
animate delete all
mol addfile $coords waitfor all

proc wrap {selText pdbOut lx ly lz} {
	set hx [expr 0.5*$lx]
	set hy [expr 0.5*$ly]
	set hz [expr 0.5*$lz]

	set sel [atomselect top $selText]
	set rAtom [$sel get {x y z}]

	set rAtomNew {}
	foreach r $rAtom {
		foreach {x y z} $r {break}

		set x [expr $x + 10*$lx]
	
		set x [expr $x-$lx*floor($x/$lx + 0.5)]
		set y [expr $y-$ly*floor($y/$ly + 0.5)]
		set z [expr $z-$lz*floor($z/$lz + 0.5)]
		lappend rAtomNew [list $x $y $z]
	}
	$sel set {x y z} $rAtomNew

	$sel writepdb $pdbOut

	$sel delete
	mol delete top
}

wrap $selText $pdbOut $lx $ly $lz
exit



