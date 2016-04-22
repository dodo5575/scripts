# Extract the coordinates of the chosen atoms.
# Use with: vmd -dispdev text -e getCoordinates.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

foreach ion {pot chl} {
    foreach nuc {ade cyt thy gua} {
	set rad 4.0
	set extraCount 20
	set selText "segname ADNA"
	#Input:
	set psf ../just_${nuc}_${ion}.psf
	set coords ../just_${nuc}_${ion}.pdb
	#Output:
	set outFile just_${nuc}_${ion}.txt

	mol load psf $psf
	mol addfile $coords
	set sel [atomselect top $selText]
	foreach zero {0} {set pos [$sel get {x y z}]}
	set n [$sel num]
	$sel delete
	mol delete top

	set out [open $outFile w]
	foreach r $pos {
	    foreach {x y z} $r {break}
	    puts $out "$x $y $z"
	    for {set i 0} {$i < $extraCount} {incr i} {
		set xi [expr {$x + (2.0*rand()-1.0)*$rad}]
		set yi [expr {$y + (2.0*rand()-1.0)*$rad}]
		set zi [expr {$z + (2.0*rand()-1.0)*$rad}]
		puts $out "$xi $yi $zi"
	    }
	}
	close $out

	puts "Wrote $n coordinates to $outFile."
    }
}
exit

