# Author: Jeff Comer <jcomer2@illinois.edu>
# Input:
set psf sio.psf
set pdb sio_ready.pdb
# Output:
set outPsf sio_ready.psf
set outMarkPdb sio_all.pdb

# Set the charges.
mol load psf $psf pdb $pdb
set all [atomselect top all]
set sil [atomselect top "type SI"]
$sil set charge 2.4
set oxy [atomselect top "type O"]
$oxy set charge -1.2
$oxy set type OSI
$all writepsf $outPsf

# Mark all atoms.
$all set beta 1.0
$all set occupancy 1.0
$all writepdb sio_all.pdb
exit
