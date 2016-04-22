set backbone "C1' H1' C2' H2' H2'' C3' O3' H3' C4' O4' H4' C5' O5' H5' H5'' O1P O2P P"
set selText "resname ADE GUA"
set psf ../distro_TA.100.psf
#set coor ../output/distro_TA.1000_run5.restart.coor
#set coor ../output/distro_TA.100_restrain1.restart.coor
set coor ../output/distro_TA.100_run5.restart.coor

#set psf ../distro_CG.100.psf
#set coor ../output/distro_CG.100_run2.restart.coor

mol load psf $psf namdbin $coor
set sel [atomselect top $selText]
set resList [lsort -unique -integer [$sel get residue]]
set resNum [llength $resList]

set resSelList {}
foreach res $resList {
    lappend resSelList [atomselect top "($selText) and residue $res"]
}

set resNum1 [expr {$resNum-1}]
for {set i 0} {$i < $resNum} {incr i } {
    set resRmsd($i) 0.0
}

for {set i 0} {$i < $resNum} {incr i} {
    set resRmsd($i) 0.0
    set posList [[lindex $resSelList $i] get {x y z}]

    for {set j 0} {$j < $resNum} {incr j} {
	set trans [measure fit [lindex $resSelList $i] [lindex $resSelList $j] weight mass]
	[lindex $resSelList $i] move $trans
	set rmsd [measure rmsd [lindex $resSelList $i] [lindex $resSelList $j] weight mass]
	set resRmsd($i) [expr {$resRmsd($i) + $rmsd}]
	[lindex $resSelList $i] set {x y z} $posList
    }
}

# for {set i 0} {$i < $resNum1} {incr i} {
#     for {set j [expr {$i+1}]} {$j < $resNum} {incr j} {
# 	set rmsd [measure rmsd [lindex $resSelList $i] [lindex $resSelList $j]]

# 	set resRmsd($i) [expr {$resRmsd($i) + $rmsd}]A
# 	set resRmsd($j) [expr {$resRmsd($j) + $rmsd}]
#     }
# }

puts ""
puts "Results"
for {set i 0} {$i < $resNum} {incr i} {
    puts "[lindex $resList $i] $resRmsd($i)"
}

foreach s $resSelList {
    $s delete
}
mol delete top
exit


