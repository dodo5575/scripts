set mole {2}
set angle 0

set pi [expr 4.0*atan(1.0)]
set a [expr $pi*$angle/180.0]
set center {0 0 0}

foreach m $mole {
    set numReps [molinfo $m get numreps]

    for {set i 0} {$i < $numReps} {incr i} {
	mol clipplane status 0 $i $m 2
	mol clipplane center 0 $i $m $center
	mol clipplane normal 0 $i $m [list [expr cos($a)] [expr sin($a)] 0]
	mol clipplane color 0 $i $m {0.75 0.75 0.75}
    }
}
