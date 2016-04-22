# mapWham.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

# Parameters:
set node 40
set tol 0.001
set waterN 1305
set springK {0.8 0.8 4.0}; # in kcal mol^1 A^-2
set cutTime 200000
set farZ0 18
set farZ1 19
# Input:
set simFile phos_pos_cen.txt
set inPrefix log/chemx_map_tclforce
set inSuffix .log.force
# Output:
set outPrefix pmf3d_tcl_high1

source whamPmf.tcl

set nx $node
set ny $node
set nz $node
set kBT 0.5862292; # in kcal/mol (at 295 K)
set waterVol 29.922; # in A^3
set avogadroNum 6.02214179e23
set volume [expr $waterN*$waterVol]

set kappa {}
foreach k $springK {
    lappend kappa [expr $k/$kBT]; # spring constant in k_B*T/A^2
}

# Read the index file.
set m [readData $simFile]
set simIndex {}
set simPos {}
set simZ {}
foreach item $m {
    lappend simIndex [lindex $item 0]
    lappend simPos [lrange $item 1 3]
}
#set simIndex [lrange $simIndex 0 3]

# Read the position data.
set data {}
set big 1e20
set x1 [expr -$big]
set x0 $big
set y1 [expr -$big]
set y0 $big
set z1 [expr -$big]
set z0 $big
set comb {}
foreach s $simIndex {
    set dataFile ${inPrefix}${s}${inSuffix}
    
    set r [readPositions $dataFile $cutTime]
    puts "Loaded [llength $r] data points from $dataFile."

    foreach p $r {
	foreach {x y z} $p {break}
	if {$x < $x0} {set x0 $x}
	if {$x > $x1} {set x1 $x}
	if {$y < $y0} {set y0 $y}
	if {$y > $y1} {set y1 $y}
	if {$z < $z0} {set z0 $z}
	if {$z > $z1} {set z1 $z}
    }
    lappend data $r
}

set nodes [list $nx $ny $nz]
set nodes [list 1 1 40]
set pmf [wham3d data $simPos $kappa $nodes 0.001]
set r [lindex $pmf 0]
set u [lindex $pmf 1]
set dx [expr ([lindex $r end 0]-[lindex $r 0 0])/($nx-1)]
set dy [expr ([lindex $r end 1]-[lindex $r 0 1])/($nx-1)]
set dz [expr ([lindex $r end 2]-[lindex $r 0 2])/($nx-1)]

set grid(nx) $nx
set grid(ny) $ny
set grid(nz) $nz
set grid(size) [expr $nx*$ny*$nz]
set grid(origin) [lindex $r 0]
set grid(delta) [list [list $dx 0 0] [list 0 $dy 0] [list 0 0 $dz]]
set grid(data) [reverseIndex u $nodes]

writeDx grid $outPrefix.dx
