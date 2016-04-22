# This script will remove all residues in the selection
# from psf and pdf files.
# Use with: vmd -dispdev text -e removeResidues.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set size 8
# Parameters:
set selText "name POT and abs(z) > 60 and x^2+y^2>15^2"
#Input:
set psf cont_dna_${size}more.psf
set pdb cont_dna_${size}more.pdb
#Output:
set finalPsf cont_dna_${size}more_neutral.psf
set finalPdb cont_dna_${size}more_neutral.pdb

set searchMax 1000
# Load the system.
mol load psf $psf pdb $pdb
set all [atomselect top all]

# Get the charge.
set q [measure sumweights $all weight charge]
set nRemove [expr int(floor($q+0.5))]

set sel [atomselect top $selText]
foreach zero {0} {set select [$sel get {segname resid name}]}
set nSelect [$sel num]
puts ""
puts "Will delete $nRemove ions from a set of $nSelect ions with ($selText)"

# Determine the total charge.
set charge [measure sumweights $all weight charge]
puts "The total charge is $charge"
$sel delete
$all delete

# Choose ions randomly from the selection.
set remove {}
set done 0
set count 0
while {$done < $nRemove && $count < $searchMax} {
    set randIndex [expr int(rand()*$nSelect)]
    set id [lindex $select $randIndex]
    
    # Don't add the same ion twice.
    if {[lsearch $remove $id] < 0} {
	lappend remove $id
	incr done
    }

    incr count
}

if {$done < $nRemove} {
    puts "Warning: Only $done ions were found."
    puts "Expand your selection or increase searchMax."
}
puts "Deleting $done ions."

package require psfgen 1.3
resetpsf

readpsf $psf
coordpdb $pdb

# Delete the selection.
foreach atom $remove {
    delatom [lindex $atom 0] [lindex $atom 1] [lindex $atom 2]
}

writepsf $finalPsf
writepdb $finalPdb
puts ""
puts "$done atoms were deleted."

# Determine the total charge of the final system.
mol delete top
mol load psf $finalPsf pdb $finalPdb
set all [atomselect top all]
set chargeNew [measure sumweights $all weight charge]
puts "Initial charge $charge"
puts "Final charge $chargeNew"
$all delete
mol delete top
exit



