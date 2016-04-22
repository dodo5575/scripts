# Get the system size using measure minmax.
# to use: vmd -dispdev text -e getSize.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set probeRadList {1 1.5 2 2.5 3 4 5}
set tol 0.2
set dx 1.0
set x0 -10
set x1 10
set y0 -10
set y1 10
set z0 8
set z1 20
set atomTextList {"type SI" "type OSI"}
set atomRadList {2.147 1.75}
set probeText "z < $z1-2"
# Input:
set name middling
set psf ${name}_surface.psf
set coord ${name}_surface_cen.pdb
# Output:
set outFile rough_${name}.dat

proc writeSurface {nodeX nodeY nodeZ outFile} {
    set out [open $outFile w]
    foreach x $nodeX y $nodeY z $nodeZ {
	puts $out "C $x $y $z"
    }
    close $out

    return
}

proc getSurfaceZ {x y z0 z1 atomTextList atomRadList probeIndex tol} {
    set za $z0
    set zb $z1
    set diff [expr $zb - $za]

    # Assume that there is a collision at z0 and none at z1.
    set probeSel [atomselect top "index $probeIndex"]
    while {$diff > $tol} {
	set zm [expr 0.5*($za+$zb)]
	$probeSel set {x y z} [list [list $x $y $zm]]

	# Look for a collision between the probe and the surface atoms.
	set collision 0
	foreach text $atomTextList rad $atomRadList {
	    set s [atomselect top "index $probeIndex and (not within $rad of (($text) and not index $probeIndex))"]
	    
	    if {[$s num] > 0} {
		set collision 1
		$s delete
		break
	    }
	    $s delete
	}

	# The old binary search.
	if {$collision} {
	    set zb $zm
	} else {
	    set za $zm
	}
	set diff [expr $zb - $za]
    }

    return [expr 0.5*($za+$zb)]
}

proc getRoughness {nodeZ} {
    set n [llength $nodeZ]

    set meanZ 0.0
    foreach z $nodeZ {
	set meanZ [expr $meanZ + $z]
    }
    set meanZ [expr $meanZ/$n]

    set roughAbs 0.0
    set roughRms 0.0
    foreach z $nodeZ {
	set roughAbs [expr $roughAbs + abs($z-$meanZ)]
	set roughRms [expr $roughRms + ($z-$meanZ)*($z-$meanZ)]
    }
    set roughAbs [expr $roughAbs/$n]
    set roughRms [expr sqrt($roughRms/$n)]
    
    return [list $meanZ $roughAbs $roughRms]
}

# Load the system.
mol load psf $psf
mol addfile $coord

# Select a probe atom.
set probeSel [atomselect top $probeText]
set probeIndex [lindex [$probeSel get index] 0]
$probeSel delete

set nx [expr int(ceil([expr ($x1-$x0)/$dx]))]
set ny [expr int(ceil([expr ($y1-$y0)/$dx]))]
puts "$nx by $ny nodes."
set nodeX {}
set nodeY {}
for {set ix 0} {$ix < $nx} {incr ix} {
    for {set iy 0} {$iy < $ny} {incr iy} {
	lappend nodeX [expr $x0 + $ix*$dx]
	lappend nodeY [expr $y0 + $iy*$dx]
    }
}


set out [open $outFile w]
foreach probeRad $probeRadList {
    # Add the probe radius to the atom radii.
    set atomRList {}
    foreach r $atomRadList {
	lappend atomRList [expr $r + $probeRad]
    }

    set nodeZ {}
    foreach x $nodeX y $nodeY {
	lappend nodeZ [getSurface $x $y $z0 $z1 $atomTextList $atomRList $probeIndex $tol]
    }

    set rough [getRoughness $nodeZ]
    puts "roughness: $probeRad $rough"
    puts $out "$probeRad [lindex $rough 2]"
}
close $out

writeSurface $nodeX $nodeY $nodeZ ${name}.xyz


mol delete top
exit
