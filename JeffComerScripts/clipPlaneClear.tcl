# Author: Jeff Comer <jcomer2@illinois.edu>
set molList [molinfo list]

foreach m $molList {
    set numReps [molinfo $m get numreps]

    for {set i 0} {$i < $numReps} {incr i} {
	mol clipplane status 0 $i $m 0
    }
}
