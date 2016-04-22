# Get the system size using measure minmax.
# to use: vmd -dispdev text -e getSize.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set selText "segname PRTA PRTB and within 4.0 of segname ADNA BDNA"
# Input:
set psf bam_spec_sys.psf
set pdb bam_spec_sys.pdb


mol load psf $psf pdb $pdb
set sel [atomselect top $selText]

puts [lsort -unique -integer [$sel get resid]]

$sel delete
mol delete top
exit



