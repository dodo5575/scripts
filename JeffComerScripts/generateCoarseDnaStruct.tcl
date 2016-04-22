# generateCoarseDnaStruct.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set pdb cg_hairpin_correct.pdb
set segName HAIR
set topoFile cg_dna.top 
set outPrefix cg_hairpin_built

package require psfgen
resetpsf
topology $topoFile

segment $segName  {
    pdb $pdb
    auto angles dihedrals
}
coordpdb $pdb $segName

writepdb $outPrefix.pdb
writepsf $outPrefix.psf

exit



