# Align restart coor for production simulation with SMD
# Usage: vmd -dispdev text -e AlignRestartCoor.tcl -args psf coorPrefix 
# Chen-Yu Li    cli56@illinois.edu
# 2014/6/5


set psf [lindex $argv 0] 
set coorPrefix [lindex $argv 1]

set inCoor "${coorPrefix}.restart.coor" 
set outCoor "${coorPrefix}_AlignRestartCoor.restart.coor"


mol load psf $psf namdbin $inCoor

set all [atomselect top all]
set dna [atomselect top nucleic]

$all moveby [vecinvert [measure center $dna]]

animate write namdbin $outCoor

exit
