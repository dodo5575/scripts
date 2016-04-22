 # This script will remove all residues in the selection
# from psf and pdf files.
# Use with: vmd -dispdev text -e removeResidues.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set angle 20.0
set diamX 3.0; # in nanometers
set diamY 2.6; # in nanometers
set parabLen 6.0; # in nanometers

#Input:
set psf neutravidin_block.psf
set pdb neutravidin_block.pdb
#Output:
set finalpsf neutravidin_pore.psf
set finalpdb neutravidin_pore.pdb

set pi [expr 4.0*atan(1.0)]
set gamma [expr $angle*$pi/180.0]
set m [expr tan($gamma)]

set a [expr {5.0*$diamX}]
set b [expr {5.0*$diamY}]
#set selText "sqrt((x/($a + abs(z)*$m))^2 + (y/($b + abs(z)*$m))^2) < 1"

set z0 [expr {5.0*$parabLen}]
set a1 [expr {$a - 0.5*$m*$z0}]
set b1 [expr {$b - 0.5*$m*$z0}]
set g [expr {0.5*$m/$z0}]

set selTextBig "sqrt((x/($a1 + abs(z)*$m))^2 + (y/($b1 + abs(z)*$m))^2) < 1"
set selTextSmall "sqrt((x/($a + $g*z^2))^2 + (y/($b + $g*z^2))^2) < 1"
set selText "(abs(z) > z0 and $selTextBig) or (abs(z) < $z0 and $selTextSmall)"
#set selText $selTextBig

# Obtain the {segid resid name} for the selection.
mol load psf $psf pdb $pdb
set sel [atomselect top $selText]
foreach null {0} {set atomList [lsort -unique [$sel get {segname resid name}]]}
set nAtoms [$sel num]
set nResidues [llength $atomList]
$sel delete

package require psfgen 1.3
resetpsf

readpsf $psf
coordpdb $pdb

# Delete the selection.
foreach atom $atomList {
    delatom [lindex $atom 0] [lindex $atom 1] [lindex $atom 2]
}

writepsf $finalpsf
writepdb $finalpdb
puts ""
puts "$nAtoms atoms were deleted."
puts "$nResidues residues were deleted."
mol delete top
exit



