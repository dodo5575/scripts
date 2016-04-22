# Author: Jeff Comer <jcomer2@illinois.edu>
source fixMutantDna1.tcl

#set tempText "segname ADNA and resid 1 to 58"
#set tempText "segname BDNA and resid 54 to 111"

set tempPdb first.pdb
set tempText "segname ADNA and resid 1 to 57"
set targPdb first.pdb
set targText "segname ADNA and resid 2 to 58"
set outPdb shiftA.pdb
fixMutantDna $tempPdb $tempText $targPdb $targText $outPdb

set tempPdb first.pdb
set tempText "segname BDNA and resid 55 to 111"
set targPdb first.pdb
set targText "segname BDNA and resid 54 to 110"
set outPdb shiftB.pdb
fixMutantDna $tempPdb $tempText $targPdb $targText $outPdb

mol load psf nw_trap2.0_0.1M.psf pdb shiftA.pdb
set sel [atomselect top "segname ADNA"]
foreach quiet {0} { set pos [$sel get {x y z}] }
mol addfile shiftB.pdb
$sel set {x y z} $pos
$sel delete

set all [atomselect top all]
$all writepdb shiftAB.pdb
$all delete
mol delete top

exit
