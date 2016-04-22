# Extract the coordinates of the chosen atoms.
# Use with: vmd -dispdev text -e getCoordinates.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set selText "segname HAIR and resid 103 to 126"
set phosphateAtoms {P O3' O5' O1P O2P}
set backboneAtoms {C1' H1' C2' H2' H2'' C3' H3' C4' H4' O4' C5' H5' H5''}
set exclude "name $phosphateAtoms $backboneAtoms"
set resnameList {ADE THY CYT GUA}

#Input:
set psf hairpin.psf
set coords hairpin.pdb
#Output:
set outFile base_atoms.txt

mol load psf $psf
mol addfile $coords

set out [open $outFile w]
foreach name $resnameList {
    set sel [atomselect top "($selText) and resname $name"]
    set res [lindex [$sel get resid] 0]
    $sel delete

    set sel [atomselect top "($selText) and not ($exclude) and resid $res"]
    set atomList [$sel get name]

    puts $out CGBEGIN
    foreach a $atomList {
	puts $out "$name\t$a\t0"
	puts "$name $a"
    }
    puts $out CGEND
}
close $out
exit



