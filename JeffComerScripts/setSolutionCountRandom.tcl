# Author: Jeff Comer <jcomer2@illinois.edu>
set sys [lindex $argv 0]
set setWater [lindex $argv 1]

set s0 15
set z0 60
#set includeText "abs(z)>$z0 and x^2+y^2>$s0^2"
set includeText "all"
set waterText "name OH2"
set ionConcList {0.1 0.1}
set ionTextList {"name POT" "name CLA"}
#Input:
set psf pore6_${sys}_conc.psf
set coor output/pore6_${sys}_eq2.restart.coor
#Output:
set finalPsf pore6_${sys}_water${setWater}.psf
set finalPdb pore6_${sys}_water${setWater}.pdb
set tempPdb tmp.pdb

# Obtain the {segid resid name} for the selection.
mol load psf $psf
mol addfile $coor

# Get the water.
set waterSel [atomselect top $waterText]
set waterN [$waterSel num]
set waterRemove [expr $waterN - $setWater]
puts "\n\nNumber of waters: $waterN"
puts "Deleting $waterRemove waters."
$waterSel delete

# Get the ions.
set ionRemoveList {}
foreach s $ionTextList c $ionConcList {
    set sel [atomselect top $s]
    set n [$sel num]
    
    set conc [expr 55.523*$n/$waterN]
    set n0 [expr $c/55.523*$waterN]
    set n1 [expr $c/55.523*$setWater]
    set remove [expr int(floor($n0 - $n1))]

    lappend ionRemoveList $remove

    puts "\nNumber of $s: $n, $n0 are free"
    puts "Concentration of $s: $conc mol/kg"
    puts "Deleting $remove ions of type $s."
    $sel delete
}

# Put the water in with the ions for elegance.
lappend ionTextList $waterText
lappend ionRemoveList $waterRemove

# Write the pdb.
set all [atomselect top all]
$all writepdb $tempPdb
set pdb $tempPdb

# Choose atoms to delete randomly.
set violators {}
foreach s $ionTextList remove $ionRemoveList {
    set sel [atomselect top "($s) and ($includeText)"]
    set nSel [$sel num]
    puts "Modifying a set of $nSel $s."
    if {$nSel < $remove} {
	puts "ERROR! Not enough $s. Reduce s0 and z0."
	exit
    }
    
    set w0 [$sel get {segname resid}]
    set w1 {}
    foreach i $w0 {
	lappend w1 [list [lindex $i 0] [lindex $i 1] [expr rand()]]
    }
    set w [lsort -real -index 2 -increasing $w1]
    set violators [concat $violators [lrange $w 1 $remove]]

    $sel delete
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
