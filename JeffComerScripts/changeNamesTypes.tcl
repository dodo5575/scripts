# This script will change the types of atoms.
# Use with: vmd -dispdev text -e changeTypes.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set changeList {{"type \"O.*\"" OSI} {"type \"SI.*\"" SI}}
set changeNameList {{"type \"O.*\"" O} {"type \"SI.*\"" SI}}
#Input:
set psf silica_surf_no_bonds.psf
set pdb silica_surf.pdb
#Output:
set finalName silica_surf_changed

mol load psf $psf pdb $pdb
set all [atomselect top all]
puts "Initial types: [lsort -unique [$all get type]]"
puts "Initial names: [lsort -unique [$all get name]]"

# Change the types.
foreach pair $changeList {
    set inText [lindex $pair 0]
    set outType [lindex $pair 1]
    set sel [atomselect top $inText]
    $sel set type $outType
    $sel delete
}

# Change the names.
foreach pair $changeNameList {
    set inText [lindex $pair 0]
    set outType [lindex $pair 1]
    set sel [atomselect top $inText]
    $sel set name $outType
    $sel set element $outType
    $sel delete
}

$all writepsf $finalName.psf
$all writepdb $finalName.pdb
puts "Final types: [lsort -unique [$all get type]]"
puts "Final names: [lsort -unique [$all get name]]"
$all delete
exit



