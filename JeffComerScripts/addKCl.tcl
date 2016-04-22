# Add ions (KCl) of given ionic strength.
# vmd -dispdev text -e addIons.tcl

# Parameters:
set conc 1.0
# Input:
set psf ssDna_sol.psf
set pdb ssDna_sol.pdb
# Output:
set prefixIons periodic_single_1M

mol load psf $psf pdb $pdb

# Compute the number of ions to add.
set posSel [atomselect top "name POT"]
set posNum [$posSel num]
set negSel [atomselect top "name CLA"]
set negNum [$negSel num]
set sel [atomselect top "name OH2"]
set nw [$sel num]
set alpha 55.523
set ni [expr int(floor($conc*$nw/($alpha + 2.0*$conc)+0.5))]
$sel delete

# Get the charge.
set all [atomselect top all]
set other [atomselect top "not name POT CLA"]
set nq [expr int(floor([measure sumweights $other weight charge]+0.5))]
set nPot [expr $ni - $posNum - $nq]
set nChl [expr $ni - $negNum]
puts "posNum0: $posNum"
puts "posNumQ: $nq"
puts "posNum1: $nPot"
puts "negNum0: $negNum"
puts "negNum1: $nChl"
$all delete
$other delete
mol delete top

package require autoionize
autoionize -psf $psf -pdb $pdb -nions [list "POT $nPot" "CLA $nChl"] -o $prefixIons -from 2.0 -between 2.0

exit
