# This script will remove all residues in the selection
# from psf and pdf files.
# Use with: vmd -dispdev text -e removeResidues.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set sys cgg
set killList {}

lappend killList [list "name POT and abs(z) > 25" 6]
lappend killList [list "nucleic" 1000]

#Input:
set psf triplet_${sys}_ions.psf
set pdb triplet_${sys}_run0-3.pdb
#Output:
set outName triplet_none_130cM

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

set finalpsf $outName.psf
set finalpdb $outName.pdb

# Obtain the {segid resid name} for the selection.
mol load psf $psf pdb $pdb

# Get the atoms to delete.
set deleteList {}
foreach kill $killList {
    foreach {selText num} $kill { break }
    set sel [atomselect top $selText]
    set selNum [$sel num]

    if {$num > $selNum} {
	puts "WARNING! Could not delete $num `$selText'. Deleting the $selNum that exist."
	set num $selNum
    }
    set last [expr {$num - 1}]
    
    set del [lrange [scramble [$sel get {segname resid name}]] 0 $last]
    $sel delete

    set deleteList [concat $deleteList $del]
}

package require psfgen 1.3
resetpsf

readpsf $psf
coordpdb $pdb

# Delete the selection.
foreach atom $deleteList {
    delatom [lindex $atom 0] [lindex $atom 1] [lindex $atom 2]
}

writepsf $finalpsf
writepdb $finalpdb
puts ""
puts "[llength $deleteList] atoms were deleted."
mol delete top
exit
