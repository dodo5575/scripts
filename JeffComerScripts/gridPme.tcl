# jcomer2@illinois.edu
source $env(HOME)/scripts/useful.tcl 
vmdargs gridPme.tcl psf pdb xsc selText pmeGridSizeX pmeGridSizeY pmeGridSizeZ outGrid

package require pmepot

set pmeSize [list $pmeGridSizeX $pmeGridSizeY $pmeGridSizeZ]

mol load psf $psf pdb $pdb
set sel [atomselect top $selText]
puts "PME atoms: [$sel num]"

pmepot -xscfile $xsc -ewaldfactor 0.25 -grid $pmeSize -sel $sel -dxfile $outGrid

$sel delete
mol delete top
exit
