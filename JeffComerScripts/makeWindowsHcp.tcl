# mapWham.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set startInd 0
set spring {0.08 0.08 0.08}
set lx 30
set ly 30
set lz 50
set poreRadius 7.0
set poreLength 26.0
set len 3
set pdb basepair_at_neg0.pdb
set excludeRad 1.5
set excludeFile basepair.txt
# Output:
set outFile win_basepair.txt
set outPointsFile win_basepair.xyz

source vector.tcl
set fy [expr sqrt(3.0)/2.0]
set fz [expr sqrt(2.0/3.0)]

set ox [expr -0.5*$lx]
set oy [expr -0.5*$ly]
set oz [expr -0.5*$lz]

set nx [expr floor($lx/$len)]
set ny [expr floor($ly/$len)]
set nz [expr floor($lz/$len)]

set posList {}
for {set iz 0} {$iz < $nz} {incr iz} {
    # Shift each layer to give hcp symmetry.
    if {$iz % 2 == 0} {
	set dxz 0
	set dyz 0
    } else {
	set dxz [expr 0.5*$len]
	set dyz [expr sqrt(3.0)/6.0*$len]
    }

    for {set iy 0} {$iy < $ny} {incr iy} {
	# Shift each row to form a hexagonal lattice.
	if {$iy % 2 == 0} {
	    set dxy 0
	} else {
	    set dxy [expr 0.5*$len]
	}
	
	for {set ix 0} {$ix < $nx} {incr ix} {
	    set x [expr $ox + $ix*$len + $dxy + $dxz]
	    set y [expr $oy + $iy*$fy*$len + $dyz]
	    set z [expr $oz + $iz*$fz*$len]
	    lappend posList [list $x $y $z]
	}
    }
}

# Find the center node.
set closePos [lindex $posList 0]
set closeDist [vecLength2 $closePos]
foreach pos $posList {
    set dist [vecLength2 $pos]
    if {$dist < $closeDist} {
	set closePos $pos
	set closeDist $dist
    }
}
puts "Shifting by $closePos."

# Shift the center node to the origin.
set cenList {}
foreach pos $posList {
    lappend cenList [vecSub $pos $closePos]
}

# Read the exclusion file.
set in [open $excludeFile r]
set excludePos {}
while {[gets $in line] >= 0} {
    if {[string length $line] < 2} { continue }
    if {[string match "\#*" $line]} { continue }
    set tok [concat $line]

    lappend excludePos $tok
}
close $in
puts "Read [llength $excludePos] exclusion sources."

# Filter the points.
set goodList {}
foreach pos $cenList {
    set s [expr sqrt([lindex $pos 0]*[lindex $pos 0] + [lindex $pos 1]*[lindex $pos 1])]

    # Enforce the geometric constraints.
    if {$s > $poreRadius || abs([lindex $pos 2]) > 0.5*$poreLength} { continue }
    # Determine whether this point is excluded.
    set excluded 0
    foreach r $excludePos {
	set d [vecLength [vecSub $pos $r]]
	if {$d < $excludeRad} { 
	    set excluded 1
	    break
	}
    }

    if {!$excluded} {
	lappend goodList $pos
    }
}
#set goodList $cenList
puts "Kept [llength $goodList] of [llength $cenList] points."

set out [open $outFile w]
set outPoints [open $outPointsFile w]
# Write the points.
puts $outPoints [llength $goodList]
puts $outPoints "window centers"
set j 0
foreach pos $goodList {
    set ind [expr $j + $startInd]
    foreach {x y z} $pos { break }
    foreach {kx ky kz} $spring { break }

    puts $out "$ind $x $y $z $kx $ky $kz"
    puts $outPoints "C $x $y $z"
    incr j
}
close $out
close $outPoints
