set pot [atomselect top "name POT"]
$pot set radius 1.76375
$pot delete

set cla [atomselect top "name CLA"]
$cla set radius 2.27
$cla delete

set si [atomselect top "resname SIO2 SIO and name \"SI.*\""]
$si set radius 2.1475
$si delete

set osi [atomselect top "resname SIO2 SIO and name \"O.*\""]
$osi set radius 1.75
$osi delete
