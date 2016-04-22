# Author: Jeff Comer <jcomer2@illinois.edu>

# Parameters:
set center [veczero]
set normal {0.0 1.0 0.0}
set selText "resname SIN"
set resolution 8
# Input:
set psf pore+dna-all.psf
set pdb conform4_run6.pdb


mol load psf $psf pdb $pdb
source solidClipPlane.tcl

writeClipSurf $center $normal top $selText $resolution pore

exit



