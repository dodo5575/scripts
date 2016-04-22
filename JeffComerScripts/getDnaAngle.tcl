# Author: Jeff Comer <jcomer2@illinois.edu>

foreach pore {18-18 20-20 23-23 24-18 22-17} {
    set frames 5
    set psf pore${pore}_ions.psf
    set prefix scaled_pore${pore}_ions

    mol load psf $psf
    set sel0 [atomselect top "(segname ADNA and resid 1) or (segname BDNA and resid 111)"]
    set sel1 [atomselect top "(segname ADNA and resid 20) or (segname BDNA and resid 92)"]

    set pi [expr {4.0*atan(1.0)}]

    for {set i 0} {$i < $frames} {incr i} {
	animate delete all
	mol addfile ${prefix}${i}.pdb

	set r0 [measure center $sel0 weight mass]
	set r1 [measure center $sel1 weight mass]
	set d [vecsub $r1 $r0]
	set dmag [veclength $d]
	set dz [lindex $d 2]
	set angle [expr {180.0/$pi*acos($dz/$dmag)}]

	puts "ANGLE: $pore $i $angle"
    }

    $sel0 delete
    $sel1 delete
    mol delete top
}
exit
