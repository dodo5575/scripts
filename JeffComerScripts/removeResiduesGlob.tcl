# This script will remove all residues in the selection
# from psf and pdf files.
# Use with: vmd -dispdev text -e removeResidues.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set pdbList [glob struct/*.pdb]
set selText "not ((segname ADNA and resid 30 31 32) or (segname BDNA and resid 29 30 31))"
#Input:
set psf trap_periodic.psf
# Output:
set outSuffix -3

proc lnmatch {lst reg} {
    set ret {}
    foreach item $lst {
	if {![regexp $reg $item]} {lappend ret $item}
    }
    return $ret
}

proc trimExtension {name} {
    set ind [string last "." $name]
    return [string range $name 0 [expr $ind-1]]
}

proc trimPath {name} {
    set ind [string last "/" $name]
    return [string range $name [expr $ind+1] end]
}

set pdbList [lnmatch $pdbList ".*-.*"]
package require psfgen 1.3

foreach pdb $pdbList {
    puts "FILE: $pdb"

    set outName "[trimExtension $pdb]${outSuffix}"
    set outPsf $outName.psf
    set outPdb $outName.pdb

    # Obtain the {segid resid name} for the selection.
    mol load psf $psf pdb $pdb
    set sel [atomselect top $selText]
    foreach zero {0} {set atomList [lsort -unique [$sel get {segname resid}]]}
    set nAtoms [$sel num]
    set nResidues [llength $atomList]
    $sel delete


    resetpsf

    readpsf $psf
    coordpdb $pdb

    # Delete the selection.
    foreach atom $atomList {
	delatom [lindex $atom 0] [lindex $atom 1]
    }

    writepsf $outPsf
    writepdb $outPdb
    puts ""
    puts "$nAtoms atoms were deleted."
    puts "$nResidues residues were deleted."
    mol delete top
}

exit
