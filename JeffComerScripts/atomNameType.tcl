# to use: vmd -dispdev text -e typeName.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set selText "segname ADNA and resid 2 60"
# Input:
set psf charmm_dna.psf
set coord charmm_dna.pdb
# Output:
set outFile name-type_charmm.txt

mol load psf $psf
mol addfile $coord
set sel [atomselect top $selText]

foreach zero {0} {set tn [lsort -unique [$sel get {resname name type charge}]]}

set out [open $outFile w]
foreach n $tn {
    puts $out "[lindex $n 0] [lindex $n 1] [lindex $n 2] [lindex $n 3]"
}
close $out

$sel delete
mol delete top
exit



