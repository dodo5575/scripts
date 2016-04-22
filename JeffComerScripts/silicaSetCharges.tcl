# Author: Jeff Comer <jcomer2@illinois.edu>
# Input:
set psf sio_ready.psf
set coor sio_anneal_parab50.restart.coor
# Output:
set outPsf silica_pore40.psf
set outPdb silica_pore40.pdb

# Set the charges.
mol load psf $psf
mol addfile $coor

set all [atomselect top all]
set sil [atomselect top "type SI"]
$sil set charge 1.0
set oxy [atomselect top "type OSI"]
$oxy set charge -0.5
$oxy set type OSI
$all writepsf $outPsf
$all writepdb $outPdb
$all delete

exit
