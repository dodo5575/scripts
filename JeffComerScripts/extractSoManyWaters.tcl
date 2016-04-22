# This script will remove all residues in the selection
# from psf and pdf files.
# Use with: vmd -dispdev text -e removeResidues.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set n 264
set waterText "name OH2"
#Input:
set psf bam_spec1.psf
set pdb bam_spec1.pdb
#Output:
set finalPsf bam_water264.psf
set finalPdb bam_water264.pdb

# Obtain the {segid resid name} for the selection.
mol load psf $psf pdb $pdb

# Get water from the top of the system.
set water [atomselect top $waterText]
foreach zero {0} {
    set w [lsort -real -index 2 -decreasing [$water get {segname resid z}]]
}
$water delete

set goodIndex {}
for {set i 0} {$i < $n} {incr i} {
    set s [atomselect top "segname [lindex $w $i 0] and resid [lindex $w $i 1]"]
    foreach ind [$s get index] {
	lappend goodIndex $ind
    }
    $s delete
}
foreach zero {0} {set selText "not (index $goodIndex)"}

# Obtain the {segid resid name} for the selection.
set sel [atomselect top $selText]
foreach zero {0} {set atomList [lsort -unique [$sel get {segname resid}]]}
set nAtoms [$sel num]
set nResidues [llength $atomList]
$sel delete

package require psfgen 1.3
resetpsf

readpsf $psf
coordpdb $pdb

# Delete the selection.
foreach atom $atomList {
	delatom [lindex $atom 0] [lindex $atom 1]
}

writepsf $finalPsf
writepdb $finalPdb
puts ""
puts "$nAtoms atoms were deleted."
puts "$nResidues residues were deleted."
mol delete top
exit



