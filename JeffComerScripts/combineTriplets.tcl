# Author: Jeff Comer <jcomer2@illinois.edu>

set inName just_
set structName struct/triplet_
set dz 6.5
set angle 50.0; # Angle base normal makes with z axis
set seg DN

package require psfgen

foreach bp(0) {a t g c} {
    foreach bp(1) {a t g c} {
	foreach bp(2) {a t g c} {
	    resetpsf

	    # Build the structure with psfgen.
	    foreach b {0 1 2} {
		if {[string equal g $bp($b)] || [string equal c $bp($b)]} {
		    readpsf ${inName}gc${b}.psf
		    coordpdb ${inName}gc${b}.pdb
		} else {
		    readpsf ${inName}at${b}.psf
		    coordpdb ${inName}at${b}.pdb
		}
	    }

	    # Write the structure.
	    set name ${bp(0)}${bp(1)}${bp(2)}
	    writepsf ${structName}${name}.psf
	    writepdb ${structName}${name}0.pdb

	    # Put the basepairs in the correct positions.
	    mol load psf ${structName}${name}.psf pdb ${structName}${name}0.pdb

	    # Rotate basepairs with pyrimidines for  z > 0.
	    foreach b {0 1 2} {
		if {[string equal t $bp($b)] || [string equal c $bp($b)]} {
		    set sel [atomselect top "segname \".${seg}${b}\""]
		    puts "num: [$sel num]"
		    $sel move [transaxis x [expr {-$angle}]]
		    $sel move [transaxis z 180.0]
		    $sel move [transaxis y 180.0]
		    $sel move [transaxis x $angle]
		    $sel delete
		}
	    }

	    # Shift the first and third.
	    set sel [atomselect top "segname \".${seg}0\""]
	    $sel moveby [list 0.0 0.0 [expr {-$dz}]]
	    $sel delete
	    set sel [atomselect top "segname \".${seg}2\""]
	    $sel moveby [list 0.0 0.0 $dz]
	    $sel delete
	    
	    set all [atomselect top all]
	    $all writepdb  ${structName}${name}.pdb
	    mol delete top
	}
    }
}

exit
