# mapWham.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

#set shortRange 1
set kT 0.5862292

if {$shortRange} {
    puts "Short Range"
    set lx 30
    set ly 30
    set lz 30
    set len 2.5
    set atomFile ${nuc}.txt
    set outPointFile nodes_${nuc}_short.xyz
    set outWindowFile win_${nuc}_short.txt
    set atomRad 7
    set minRad 2.2
} else {
    puts "Long Range"
    set lx 80
    set ly 80
    set lz 80
    set len 5
    set atomFile ${nuc}.txt
    set outPointFile nodes_${nuc}_long.xyz
    set outWindowFile win_${nuc}_long.txt
    set atomRad 16.0
    set minRad 4
}

set spring [expr {$kT/($len*$len)}]

source $env(HOME)/scripts/vector.tcl
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
set in [open $atomFile r]
set atomPos {}
while {[gets $in line] >= 0} {
    if {[string length $line] < 2} { continue }
    if {[string match "\#*" $line]} { continue }
    set tok [concat $line]

    lappend atomPos $tok
}
close $in
puts "Read [llength $atomPos] exclusion sources."

set removed 0

# Filter the points.
set goodList {}
foreach pos $cenList {
    # Determine whether this point is excluded.
    set excluded 1
    foreach r $atomPos {
	set d [vecLength [vecSub $pos $r]]
	if {$d < $atomRad} { 
	    set excluded 0
	    break
	}
    }

    if {!$excluded} {
	foreach r $atomPos {
	    set d [vecLength [vecSub $pos $r]]
	    if {$d < $minRad} { 
		set excluded 1
		incr removed
		break
	    }
	}
    }

    if {!$excluded} {
	lappend goodList $pos
    }
}
#set goodList $cenList
puts "Kept [llength $goodList] of [llength $cenList] points."

set outPoint [open $outPointFile w]
set outWindow [open $outWindowFile w]
# Write the points.
puts $outPoint [llength $goodList]
puts $outPoint "window centers"
set j 0
foreach pos $goodList {
    set ind [expr $j + $start]
    foreach {x y z} $pos { break }

    puts $outPoint "C $x $y $z"
    puts $outWindow "$ind $x $y $z $spring $spring $spring"
    incr j
}
close $outPoint
close $outWindow

puts "number too close: $removed"
set start [expr $j+$start]




