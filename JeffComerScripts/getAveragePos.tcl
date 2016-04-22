# Author: Jeff Comer <jcomer2@illinois.edu>
source $env(HOME)/scripts/vmdargs.tcl

vmdargslist getAveragePos.tcl psf pdb stride selText outPdb @dcdList

puts "DCDLIST: $dcdList"

mol load psf $psf pdb $pdb
set all [atomselect top all]
set sel [atomselect top $selText]
set posList0 [$all get {x y z}]

foreach dcd $dcdList {
    mol addfile $dcd step $stride waitfor all
}

set posList [measure avpos $sel]
$all set {x y z} $posList0
$sel set {x y z} $posList
$all writepdb $outPdb

$sel delete
$all delete
mol delete top
exit
