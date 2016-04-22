# Move the center of mass to the origin.
# use with: vmd -dispdev text -e center.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>
set moveText all
set dirText0 "(segname ADNA and resid 20) or (segname BDNA and resid 43)"
set dirText1 "(segname ADNA and resid 1) or (segname BDNA and resid 62)"
set alignDir {0 0 -1}

# Input:
set psf bam_complex4.psf
set coor bam_complex4.pdb
# Output:
set finalPdb bam_complex4_align.pdb

mol load psf $psf
mol addfile $coor waitfor all

# Get the vector.
set dirSel0 [atomselect top $dirText0]
set cen0 [measure center $dirSel0 weight mass]
$dirSel0 delete

set dirSel1 [atomselect top $dirText1]
set cen1 [measure center $dirSel1 weight mass]
$dirSel1 delete

# Move stuff.
set moveSel [atomselect top $moveText]
set v [vecsub $cen1 $cen0]
set m0 [transvecinv $v]
set m1 [transvec $alignDir]
set m [transmult $m1 $m0]
$moveSel move $m
$moveSel delete

# Write the papers.
set allSel [atomselect top all]
$allSel writepdb $finalPdb
$allSel delete

mol delete top
exit



