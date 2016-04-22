source periodicCell.tcl
set a [expr {1.250483*19}]
set b [expr {1.241057*19}]
set c [expr {0.247888*116}]

set numFrames [molinfo top get numframes]
for {set i 0} {$i < $numFrames} {incr i} {
    molinfo top set frame $i
    setCell $a $b $c 90 90 90
}
