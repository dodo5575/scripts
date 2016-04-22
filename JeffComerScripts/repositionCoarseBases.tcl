# Extract the coordinates of the chosen atoms.
# Use with: vmd -dispdev text -e getCoordinates.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set posText "(resname ADE GUA and name N1) or (resname CYT THY and name N3)"
#Input:
set atomPsf hairpin.psf
set atomCoords hairpin.pdb
set coarsePdb cg_hairpin.pdb
#Output:
set outFile cg_hairpin_correct.pdb

# Get the atomistic coordinates.
mol load psf $atomPsf 
mol addfile $atomCoords
set sel [atomselect top $posText]
set res [$sel get resid]
set pos [$sel get {x y z}]
$sel delete
mol delete top

# Load the coarse grain structure.
mol load pdb $coarsePdb
set all [atomselect top all]

# Set the coordinates.
foreach i $res r $pos {
    set sel [atomselect top "resid $i and name \"\[ATCG\]B\""]
    $sel set {x y z} [list $r]
    $sel delete
}

# Write the results.
$all writepdb $outFile
$all delete

exit



