# This script will remove waters starting at the bottom
# of the system (up to z0) to obtain the nWater waters.
# Use with: vmd -dispdev text -e setWaterCount.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set s0 15
set z0 55
#set nWater 73914
set nWater 72000
set waterText "name OH2"
#Input:
set psf extra_water7226.psf
#set pdb extra_water7226.pdb
set coor output/dna_extra7226_eq1.restart.coor
#Output:
set finalPsf dna_water${nWater}.psf
set finalPdb dna_water${nWater}.pdb
set tempPdb tmp.pdb

# Obtain the {segid resid name} for the selection.
mol load psf $psf
mol addfile $coor

set all [atomselect top all]
$all writepdb $tempPdb
set pdb $tempPdb

# Find the number of waters to remove.
set water [atomselect top $waterText]
set nWater0 [$water num]
set n [expr $nWater0 - $nWater]
puts "Removing $n waters..."
$water delete

set sel [atomselect top "$waterText and abs(z)>$z0 and x^2+y^2>$s0^2"]
set nSel [$sel num]
puts "Modifying a set of $nSel waters."
if {$nSel < $n} {
    puts "ERROR! Not enough water in original system."
    exit
}

# Remove low water first.
puts "Sorting $nSel waters..."
# Sort by something random, like the fraction part of the z position.
foreach zero {0} {
    #set w [lsort -real -index 2 -increasing [$sel get {segname resid z}]]
    set w0 [$sel get {segname resid z}]
    set w1 {}
    foreach i $w0 {
	set z [expr 100*[lindex $i 2]]
	lappend w1 [list [lindex $i 0] [lindex $i 1] [expr $z-floor($z)]]
    }
    set w [lsort -real -index 2 -increasing $w1]
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

