# Add harmonic constraints to silicon nitride.
# to use: vmd -dispdev text -e constrainSilicon.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

# Parameters:
# Spring constant in kcal/(mol A^2)
set selText "(resname SIN and name \"N.*\" and within 5.0 of not resname SIN) or (resname SIN and name \"SI.*\" and within 5.049668 of not resname SIN)"
set distribText "resname SIN and not ($selText)"
set factor 2.0
# Input:
set psf pore+dna_E4V.psf
set pdb steady_4V2.pdb
# Output:
set psfOut charged.psf

mol load psf $psf pdb $pdb
set all [atomselect top all]
set totalCharge [measure sumweights $all weight charge]
puts "Total charge: $totalCharge"

set sel [atomselect top $selText]
set distrib [atomselect top $distribText]
set charge0 [measure sumweights $sel weight charge]
puts "Modifying the charge of [$sel num] atoms: $selText"
puts "Distributing to [$distrib num] atoms: $distribText"
puts "Charge of selection: $charge0"

set qSel {}
foreach q [$sel get charge] {
    lappend qSel [expr $q*$factor]
}
$sel set charge $qSel
set charge [measure sumweights $sel weight charge]
set chargeDistrib [expr -($charge-$charge0)]

set nDistrib [$distrib num]
set dq [expr $chargeDistrib/$nDistrib]

puts "Distributing charges of $dq"
set qDistrib {}
foreach q [$distrib get charge] {
    lappend qDistrib [expr $q + $dq]
}
$distrib set charge $qDistrib

set totalCharge [measure sumweights $all weight charge]
puts "Total charge now: $totalCharge"

puts "Neutralizing the whole system"
foreach x {0} {lset qDistrib 0 [expr [lindex $qDistrib 0] - $totalCharge]}
$distrib set charge $qDistrib

set totalCharge [measure sumweights $all weight charge]
puts "Total charge now: $totalCharge"

$all writepsf $psfOut
exit



