# Add harmonic constraints to silicon nitride.
# to use: vmd -dispdev text -e langevinWater.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

#foreach sys {q0 q1 q-1 q2c} {}
foreach sys {q2r q2s} {
# Parameters:
# Langevin constant in ps^-1.
set betaList {1}
set selText "resname SIO2"
# Input:
set psf anneal_dopc_${sys}.psf
set pdb anneal_dopc_${sys}.pdb
# Output:
set restFilePrefix layer_langevin_${sys}

mol load psf $psf pdb $pdb
set selAll [atomselect top all]

# Set the spring constants to zero for all atoms.
$selAll set occupancy 0.0
$selAll set beta 0.0

# Select the atoms to which langevin is applied.
set sel [atomselect top $selText]

foreach beta $betaList {
	# Set the Langevin coupling.
	$sel set beta $beta
	# Write the constraint file.
	$selAll writepdb $restFilePrefix.pdb
}

puts "\nSet langevin damping to $betaList ps^-1 for [$sel num] atoms.\n"

$sel delete
$selAll delete
mol delete top
}
exit
