# Make a pore geometry file cylindrical pore with rounded corners.
# Use with: tclsh phantomGeometry.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

# Parameters:
set dz 0.5
set cornerRadius 20.0
set poreLength 110.0
set poreRadius 10.0
# Output:
set outFile pore${poreRadius}-${poreLength}.geo

set pi [expr 4.0*atan(1.0)]
set n [expr int(ceil($poreLength/$dz))]

# Return the radius as a function of z.
proc radius {z l rp rc} {
    set hl [expr 0.5*$l]
    if {$z > $hl - $rc} {
	return [expr $rp + $rc - sqrt($rc*$rc - ($z-$hl+$rc)*($z-$hl+$rc))]
    } elseif {$z < -$hl + $rc} {
	return [expr $rp + $rc - sqrt($rc*$rc - ($z+$hl-$rc)*($z+$hl-$rc))]
    } else {
	return $rp
    }
}

# Write the geometry.
set out [open $outFile w]
for {set i 0} {$i <= $n} {incr i} {
    set z [expr $poreLength*(1.0*$i/$n - 0.5)]
    set s [radius $z $poreLength $poreRadius $cornerRadius]
    puts $out "$z $s"
}
close $out
puts "$n radius values were written."
exit



