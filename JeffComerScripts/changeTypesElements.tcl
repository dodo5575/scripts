# This script will change the types of atoms.
# Use with: vmd -dispdev text -e changeTypes.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set changeTypeList {{"type \"O.*\"" OSI} {"type \"SI.*\"" SI}}
set changeElementList {{"type \"O.*\"" O} {"type \"SI.*\"" Si}}
#Input:
set psf silica_surf_no_bonds.psf
set pdb silica_surf.pdb
#Output:
set finalName silica_surf_changed

mol load psf $psf pdb $pdb
set all [atomselect top all]
puts "Initial types: [lsort -unique [$all get type]]"
puts "Initial elements: [lsort -unique [$all get element]]"

# Change the types.
foreach pair $changeTypeList {
    set inText [lindex $pair 0]
    set outType [lindex $pair 1]
    set sel [atomselect top $inText]
    $sel set type $outType
    $sel delete
}

# Change the elements.
foreach pair $changeElementList {
    set inText [lindex $pair 0]
    set outType [lindex $pair 1]
    set sel [atomselect top $inText]
    $sel set element $outType
    $sel delete
}

$all writepsf $finalName.psf
$all writepdb $finalName.pdb
puts "Final types: [lsort -unique [$all get type]]"
puts "Final elements: [lsort -unique [$all get element]]"
$all delete
exit



