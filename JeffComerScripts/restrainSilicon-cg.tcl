# Add restraints to coarse-grained silicon.
# to use: vmd -dispdev text -e restrainSilicon-cg.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

# Input:
set psf system-ions.psf
set pdb system-ions.pdb
# Output:
# for fixing silicon
set beta 10.0
set restFile siliconRest_${beta}.pdb
# Distance between silicon beads.
set dist 5.0

#------------ silicon constraints ----#
mol load psf $psf pdb $pdb
set selAll [atomselect top all]
$selAll set occupancy 0.0
$selAll set beta 0.0

set selSIN [atomselect top "resname SIN SIO SOH"]
$selSIN set beta $beta

$selAll writepdb $restFile
exit



