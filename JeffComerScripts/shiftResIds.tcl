# Shift the resIds to start with 1.
# Use with: vmd -dispdev text -e shiftResIds.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set selText "segname ADNA"
#Input:
set psf dsDna62.psf
set pdb dsDna62.pdb
#Output:
set finalPsf dsDna62_res.psf
set finalPdb dsDna62_res.pdb

# Obtain the {segid resid name} for the selection.
mol load psf $psf pdb $pdb
set sel [atomselect top $selText]

# Find the resIds.
foreach zero {0} {set residue [lsort -unique -integer [$sel get resid]]}
set res0 [lindex $residue 0]
puts "Root residue: $res0"
$sel delete

# Select the residues.
set selList {}
foreach res $residue {
    lappend selList [atomselect top "($selText) and resid $res"]    
}

# Shift the resIds.
foreach s $selList res $residue {
    $s set resid [expr $res - $res0 + 1]
    $s delete
}

# Save the results.
set all [atomselect top all]
$all writepsf tmp.psf
$all writepdb $finalPdb

# Copy the atom records so that angles, dihedrals, etc. are not lost.
#source insertPsfAtomRecords.tcl
puts "tclsh insertPsfAtomRecords.tcl tmp.psf $psf $finalPsf"

mol delete top
exit



