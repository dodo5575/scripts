# Author: jcomer2@illinois.edu
# Parameterps:
source $env(HOME)/scripts/useful.tcl

vmdargs topographyMap.tcl psf pdb reverse outputPrefix

set startPos {0 0}
set endPos {25 25}
set startZ 18
set endZ 0
set stepZ -0.05
set len 0.25
set typeList {SI OSI}
set radiusList {2.1475 1.75}

if {$reverse > 0} {
    set stepZ [expr {-$stepZ}]
    set startZ [expr {-$startZ}]
}

set x0 [lindex $startPos 0]
set y0 [lindex $startPos 1]
set lx [expr {[lindex $endPos 0]-$x0}]
set ly [expr {[lindex $endPos 1]-$y0}]
set nx [expr {int(ceil($lx/$len))}]
set ny [expr {int(ceil($ly/$len))}]


proc probe {x y startZ endZ type rad} {
    global stepZ

    set sel [atomselect top water]
    set probeInd [lindex [$sel get index] 0]
    $sel delete

    set sel [atomselect top "index $probeInd"]
    
    for {set z $startZ} {($stepZ > 0.0 && $z < $endZ) || ($stepZ < 0.0 && $z > $endZ)} {set z [expr {$z+$stepZ}]} {
	$sel set {x y z} [list [list $x $y $z]]

	set s [atomselect top "index $probeInd and within $rad of type $type"]
	set num [$s num]
	$s delete

	if {$num > 0} {
	    # We found the surface.
	    # Stop here.
	    #puts $z
	    break
	}
    }
    $sel delete

    return $z
}

# Scan the surface.
set pdb $pdb

mol load pdb $pdb

set sel [atomselect top water]
set atomInd [lindex [$sel get index] 0]

set outX [open $outputPrefix.dat.rx w]
set outY [open $outputPrefix.dat.ry w]
set outZ [open $outputPrefix.dat.rz w]

set total [expr {$nx*$ny}]
set complete -100

for {set ix 0} {$ix < $nx} {incr ix} {
    for {set iy 0} {$iy < $ny} {incr iy} {
	set x [expr {$x0 + $lx*$ix/double($nx)}]
	set y [expr {$y0 + $ly*$iy/double($ny)}]

	# Find the closest surface point.
	set surfZ $endZ
	foreach type $typeList rad $radiusList {
	    set z  [probe $x $y $startZ $endZ $type $rad]
	    if {($stepZ < 0.0 && $z > $surfZ) || ($stepZ > 0.0 && $z < $surfZ)} {
		set surfZ $z
	    }
	}

	# Write the data.
	if {$iy == 0} {
	    puts -nonewline $outX "$x"
	    puts -nonewline $outY "$y"
	    puts -nonewline $outZ "$surfZ"
	} else {
	    puts -nonewline $outX " $x"
	    puts -nonewline $outY " $y"
	    puts -nonewline $outZ " $surfZ"
	}
    }
    puts $outX ""
    puts $outY ""
    puts $outZ ""

    # Status.
    set comp [expr {(100*$ix*$iy)/$total}]
    if {abs($complete - $comp) >= 5} {
	puts "$comp percent complete"
	set complete $comp
    }
}

close $outX
close $outY
close $outZ

exit
