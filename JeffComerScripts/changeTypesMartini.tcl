# This script will change the types of atoms.
# Use with: vmd -dispdev text -e changeTypes.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>
set conc 1.0

set changeList {{"type P" P4}}
set posIon SOD
set negIon CL
#Input:
set psf cg_ions_${conc}M.psf
set pdb cg_ions_${conc}M.pdb
#Output:
set finalPsf cg_ions_${conc}M_martini.psf

mol load psf $psf pdb $pdb
set all [atomselect top all]
puts "Initial types: [lsort -unique [$all get type]]"

set sel [atomselect top "name $posIon"]
$sel set charge 1.0
$sel delete

set sel [atomselect top "name $negIon"]
$sel set charge -1.0
$sel delete
puts "Shifted the charges of $posIon and $negIon to +1 and -1."

# Change the types.
foreach pair $changeList {
    set inText [lindex $pair 0]
    set outType [lindex $pair 1]
    set sel [atomselect top $inText]
    $sel set type $outType
    $sel delete
}

$all writepsf $finalPsf
puts "Final types: [lsort -unique [$all get type]]"
$all delete
exit



