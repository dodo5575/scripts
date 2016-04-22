# This script combines two pdb/psf pairs.
# Use with: vmd -dispdev text -e combine.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

# Input:
set psfList {ADNA_charmm.psf BDNA_charmm1.psf other.psf}
set pdbList {ADNA.pdb BDNA_charmm1.pdb other.pdb}
# Output:
set finalPsf cont_dna_charmm.psf
set finalPdb cont_dna_charmm.pdb

# Load the topology and coordinates.
package require psfgen
resetpsf
foreach psf $psfList pdb $pdbList {
    readpsf $psf
    coordpdb $pdb
}

# Write the combination.
writepdb $finalPdb
writepsf $finalPsf
exit


