set angle -90.0
set d 38
set d1 35

set pi [expr 4.0*atan(1.0)]
set t [expr $pi/6.0]
set a [expr $pi*$angle/180.0]

mol clipplane status 0 0 1 2
mol clipplane center 0 0 1 {0 0 0}
mol clipplane normal 0 0 1 [list [expr -cos($a)] [expr -sin($a)] 0]
mol clipplane color 0 0 1 {0.75 0.75 0.75}

mol clipplane status 1 0 1 2
mol clipplane center 1 0 1 [list [expr -$d*cos($a)] [expr -$d*sin($a)] 0]
mol clipplane normal 1 0 1 [list [expr cos($a+$t)] [expr sin($a+$t)] 0]
mol clipplane color 1 0 1 {0.75 0.75 0.75}

mol clipplane status 2 0 1 2
mol clipplane center 2 0 1 [list [expr -$d*cos($a)] [expr -$d*sin($a)] 0]
mol clipplane normal 2 0 1 [list [expr cos($a-$t)] [expr sin($a-$t)] 0]
mol clipplane color 2 0 1 {0.75 0.75 0.75}

mol clipplane status 3 0 1 2
mol clipplane center 3 0 1 [list [expr -$d1*sin($a)] [expr $d1*cos($a)] 0]
mol clipplane normal 3 0 1 [list [expr sin($a)] [expr -cos($a)] 0]
mol clipplane color 3 0 1 {0.75 0.75 0.75}

mol clipplane status 4 0 1 2
mol clipplane center 4 0 1 [list [expr $d1*sin($a)] [expr -$d1*cos($a)] 0]
mol clipplane normal 4 0 1 [list [expr -sin($a)] [expr cos($a)] 0]
mol clipplane color 4 0 1 {0.75 0.75 0.75}



