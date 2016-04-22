set selText "nucleic"
set beta 1.0
# Input:
set psf one_turn_at_pot.psf
set refPdb one_turn_at_pot.pdb
set fitPdb scaled_one_turn_at_pot.pdb
# Output:
set outName good

set refMol [mol load psf $psf pdb $refPdb]
set fitMol [mol load psf $psf pdb $fitPdb]

set refSel [atomselect $refMol $selText]
set fitSel [atomselect $fitMol $selText]

# Set the coordinates.
set pos [$refSel get {x y z}]
$fitSel set {x y z} $pos
$fitSel set beta $beta
puts "Set [$fitSel num] atoms."

$refSel delete
$fitSel delete

set all [atomselect $fitMol all]
$all writepsf $outName.psf
$all writepdb $outName.pdb
$all delete

mol delete $fitMol
mol delete $refMol

exit
