# Author: Jeff Comer <jcomer2@illinois.edu>
source $env(HOME)/scripts/vmdargs.tcl

vmdargslist getCoordTypeRadius.tcl psf coordFile selText typeRadiusFile outFile

# Just read a space delimited data file.
proc readData {fileName} {
    set in [open $fileName r]
    
    set r {}
    while {[gets $in line] >= 0} {
	if {[string match "#*" $line]} {continue}
	if {[string length $line] < 2} {continue}

	set tok [concat $line]
	lappend r $tok
    }

    close $in
    return $r
}

proc quiet {} {}

mol load psf $psf
mol addfile $coordFile

# Set the radii of the atoms.
set data [readData $typeRadiusFile]; quiet
foreach item $data {
    set s [atomselect top "type [lindex $item 0]"]
    if {[$s num] > 0} {
	$s set radius [lindex $item 1]
    }
    $s delete
}

# Write the positions and radii.
set sel [atomselect top $selText]
foreach zero {0} {set pos [$sel get {radius x y z}]}
set n [$sel num]
$sel delete
mol delete top

set out [open $outFile w]
foreach r $pos {
    foreach {x y z rad} $r {break}
    puts $out "$x $y $z $rad"
}
close $out

puts "Wrote $n coordinates to $outFile."
exit

