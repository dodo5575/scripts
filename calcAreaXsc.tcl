# calculate the area in each xsc and return the one with the closest to mean value
# Usage: tclsh calcAreaXsc.tcl xscName NumOfXsc output dim meanArea 
# Chen-Yu Li    cli56@illinois.edu
# 2014/6/5


source /home/cli56/scripts/Procs.tcl

set xscName [lindex $argv 0]
set NumOfXsc [lindex $argv 1]
set output [lindex $argv 2]
set out [open "${output}_calcAreaXsc.dat" w]
set dim [lindex $argv 3]
set meanArea [lindex $argv 4]

set dArea 99999

for {set i 1} {$i <= $NumOfXsc} {incr i} {

    set xsc [readExtendedSystem ${xscName}.${i}.restart.xsc]
    set x [lindex $xsc 1]
    set y [lindex $xsc 5]
    set z [lindex $xsc 9]

    switch -- $dim {
    
        "x" {
            set area [expr $y * $z]
        }
        "y" {
            set area [expr $x * $z]
        }
        "z" {
            set area [expr $x * $y]
        }
        default {
        }
    }

    puts $out [format "${xscName}.${i}.restart.xsc\t%.3f\t%.3f\t%.3f\t%.3f" $area $x $y $z]

    set dArea_tmp [expr abs($meanArea - $area)]
    if {$dArea_tmp <= $dArea} {
        set dArea $dArea_tmp
        set closeXsc "${xscName}.${i}.restart.xsc"
        set closeArea $area
        set closeX $x
        set closeY $y
        set closeZ $z
    }

}
puts $out [format "\nInput Mean Area = %.3f" $meanArea]
puts $out [format "$closeXsc\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f" $dArea $closeArea $closeX $closeY $closeZ]

close $out

puts [format "\nInput Mean Area = %.3f" $meanArea]
puts [format "$closeXsc\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f" $dArea $closeArea $closeX $closeY $closeZ]


