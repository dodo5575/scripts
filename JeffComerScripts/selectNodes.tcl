mol load psf comb.psf pdb comb.pdb

set x0 -4
set x1 2
set y0 -1
set y1 4
set z0 -5
set z1 7

set selText "segname NULL and x > $x0 and x < $x1 and y > $y0 and y < $y1 and z > $z0 and z < $z1"
set sel [atomselect top "$selText"]

foreach quiet {0} {set nodeList [$sel get resid]}
puts $nodeList
exit
