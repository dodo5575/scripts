# This script will remove all residues in the selection
# from psf and pdf files.
# Use with: vmd -dispdev text -e removeResidues.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set killList {}

lappend killList [list "name OH2" 18625]
set z0 0.0
#Input:
set psf bilayer_dopc_1M.psf
set pdb run2_dopc_wat.pdb
#Output:
set outName dopc_short

proc scramble {inList} {
    set sortList {}
    foreach item $inList {
	lappend sortList [list [expr {rand()}] $item]
    }

    # Sort based on the random number.
    set sortList [lsort -real -index 0 $sortList]

    set outList {}
    foreach item $sortList {
	lappend outList [lindex $item 1]
    }
    return $outList
}

proc sortAbsZ {inList z0} {
     set sortList {}
    foreach item $inList {
	foreach {seg res z} $item { break }
	lappend sortList [list $seg $res [expr {abs($z-$z0)}]]
    }

    # Sort based on abs(z).
    set sortList [lsort -real -decreasing -index 2 $sortList]

    return $sortList
}

set finalpsf $outName.psf
set finalpdb $outName.pdb

# Obtain the {segid resid name} for the selection.
mol load psf $psf pdb $pdb

# Get the atoms to delete.
set deleteList {}
foreach kill $killList {
    foreach {selText num} $kill { break }
    set sel [atomselect top $selText]
    set allList [sortAbsZ [$sel get {segname resid z}] $z0]
    set del [lrange $allList 0 [expr {$num-1}]]
    $sel delete

    set deleteList [concat $deleteList $del]
}

package require psfgen 1.3
resetpsf

readpsf $psf
coordpdb $pdb

# Delete the selection.
foreach atom $deleteList {
    delatom [lindex $atom 0] [lindex $atom 1]
}

writepsf $finalpsf
writepdb $finalpdb
puts ""
puts "[llength $deleteList] residues were deleted."
mol delete top
exit
