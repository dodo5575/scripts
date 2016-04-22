# Author: Jeff Comer <jcomer2@illinois.edu>
set outName cyl2a.dat
set gridFile ${outName}.grid
set lenZ 72.0
set lenX 40.0
set delZ 2.0
set delX 2.0

set lz [expr {0.5*$lenZ}]
set ls [expr {0.5*$lenX}]
set z0 0.0
set s0 0.0

set nz [expr {int(ceil($lz/$delZ))}]
set ns [expr {int(ceil($ls/$delX))}]
set ds [expr {$ls/$ns}]
set dz [expr {$lz/$nz}]

set out [open $outName w]
puts $out "$s0 $z0"
puts $out "$ds $dz"
puts $out "$ns $nz"
close $out

set out [open $gridFile w]
set j 0
for {set is 0} {$is < $ns} {incr is} {
    for {set iz 0} {$iz < $nz} {incr iz} {
	set s [expr {$s0 + $ds*$is}]
	set z [expr {$z0 + $dz*$iz}]
	puts $out "$j $s $z"
	incr j
    }
}
close $out
