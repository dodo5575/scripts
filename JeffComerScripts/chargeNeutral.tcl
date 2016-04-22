# This script will remove all residues in the selection
# from psf and pdf files.
# Use with: vmd -dispdev text -e removeResidues.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set selText0 "type SI"
set selText "type O"

#Input:
set psf  silica_shell_good0.psf
set pdb  silica_shell_good.pdb
#Output:
set finalpsf silica_shell_good.psf

# Obtain the {segid resid name} for the selection.
mol load psf $psf pdb $pdb
set all [atomselect top all]
set sel0 [atomselect top $selText0]
set sel [atomselect top $selText]

set num0 [$sel0 num]
set charge0 [measure sumweights $sel0 weight charge]
set num [$sel num]
set charge [measure sumweights $sel weight charge]

set qi [lindex [$sel get charge] 0]
set q [expr -$charge0/$num]
set err [expr ($q-$qi)/$qi]
puts "Setting the charge of `$selText' from $qi to $q"
puts "Relative error: $err"
$sel set charge $q

$all writepsf $finalpsf

$sel0 delete
$sel delete
$all delete
mol delete top
exit



