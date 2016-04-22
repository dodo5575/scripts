# This script will remove waters starting at the bottom
# of the system (up to z0) to obtain the nWater waters.
# Use with: vmd -dispdev text -e setWaterCount.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set z0 -55
#set nWater 73914
set nWater 72000
set waterText "name OH2"
#Input:
set psf extra_water7226.psf
set pdb extra_water7226.pdb
#Output:
set finalPsf dna_water${nWater}.psf
set finalPdb dna_water${nWater}.pdb

# Obtain the {segid resid name} for the selection.
mol load psf $psf pdb $pdb

# Find the number of waters to remove.
set water [atomselect top $waterText]
set nWater0 [$water num]
set n [expr $nWater0 - $nWater]
puts "Removing $n waters..."
$water delete

set sel [atomselect top "$waterText and z < $z0"]
set nSel [$sel num]
puts "Modifying a set of $nSel waters."
if {$nSel < $n} {
    puts "ERROR! Not enough water in original system."
    exit
}

# Remove low water first.
puts "Sorting $nSel waters..."
foreach zero {0} {
    set w [lsort -real -index 2 -increasing [$sel get {segname resid z}]]
}
$sel delete

puts "Collecting violators..."
# Collect the indices of the water to delete.
set violators {}
for {set i 0} {$i < $n} {incr i} {
    lappend violators [lrange [lindex $w $i] 0 1]
    #puts [lrange [lindex $w $i] 0 1]
}

package require psfgen 1.3
resetpsf

readpsf $psf
coordpdb $pdb

# Delete the selection.
foreach atom $violators {
	delatom [lindex $atom 0] [lindex $atom 1]
}

writepsf $finalPsf
writepdb $finalPdb
puts ""
puts "[llength $violators] residues were deleted."
mol delete top
exit



