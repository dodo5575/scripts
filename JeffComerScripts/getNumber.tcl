# Get the system size using measure minmax.
# to use: vmd -dispdev text -e getSize.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set selTextList {"name POT" "name CLA" "name OH2"}
# Input:
set psf bam_solv1.psf
set pdb bam_solv1.pdb


mol load psf $psf pdb $pdb
foreach s $selTextList {
    set sel [atomselect top $s]
    puts "Number of $s: [$sel num]"
    $sel delete
}
mol delete top
exit



