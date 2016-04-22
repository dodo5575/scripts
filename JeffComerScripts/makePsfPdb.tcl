# Author: Jeff Comer <jcomer2@illinois.edu>

# Input
set inName diff_grid.xyz

proc trimExtension {name} {
    set ind [string last "." $name]
    return [string range $name 0 [expr $ind-1]]
}

set prefix [trimExtension $inName]

mol new $inName
set all [atomselect top all]
$all set resid [$all get index]
$all set resname C
$all set type C
$all set name C
$all set segname NULL 
$all writepsf $prefix.psf
$all writepdb $prefix.pdb
$all delete

mol delete top
exit
