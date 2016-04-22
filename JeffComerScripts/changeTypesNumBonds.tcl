# This script will change the types of atoms.
# Use with: vmd -dispdev text -e changeTypes.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set changeList {{"type N" N} {"type SI" SI}}
#Input:
set psf trap20_bigger.psf
set pdb trap20_bigger.pdb
#Output:
set finalPsf trap20_bigger_types.psf

mol load psf $psf pdb $pdb
set all [atomselect top all]
puts "Initial types: [lsort -unique [$all get type]]"

# Change the types.
foreach pair $changeList {
    set inText [lindex $pair 0]
    set outType [lindex $pair 1]
    set sel [atomselect top $inText]
    set numList [$sel get numbonds]
    set typeList0 [$sel get type]

    set typeList {}
    foreach n $numList t $typeList0 {
	lappend typeList ${t}_${n}
    }
    $sel set type $typeList
    
    $sel delete
}

$all writepsf $finalPsf
puts "Final types: [lsort -unique [$all get type]]"
$all delete
exit
