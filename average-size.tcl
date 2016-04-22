set cutTime 10000
set in [open output/pbcPlate-nvt-langevin.xst r]
set xysum 0.0
set zsum 0.0
set count 0
while { [gets $in line] >= 0 } {
    if { [string match "#*" $line] } { continue }
    set tok [concat $line]
    if { [lindex $tok 0] < $cutTime } { continue }
    set xysum [expr {$xysum + [lindex $tok 1]}]
    set zsum [expr {$zsum + [lindex $tok 9]}]
    incr count
}
set xymean [expr {$xysum/$count}]
set zmean [expr {$zsum/$count}]
puts "xymean = $xymean"
puts "zmean = $zmean"
