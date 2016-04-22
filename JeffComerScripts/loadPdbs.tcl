# Author: Jeff Comer <jcomer2@illinois.edu>
set stride 1
set startFrame 2
set nFrames 33
# Input:
set prefix short_flow
set inDir frames

proc viewMol {molId} {
    set molList [molinfo list]
    
    foreach m $molList {
	molinfo $m set displayed off
    }
    molinfo $molId set displayed on
}

for {set i $startFrame} {$i < $nFrames} {incr i $stride} {
    mol new $inDir/$prefix.$i.pdb
    molinfo top set displayed off
}
