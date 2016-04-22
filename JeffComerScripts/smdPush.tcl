# Add harmonic constraints to silicon nitride.
# to use: vmd -dispdev text -e constrainSilicon.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

#set water 110413
#set sys s-dna
set water 91924
set sys b-dna
# Parameters:
# Spring constant in kcal/(mol A^2)
set occupancy 1.0
set selText "segname ADNA BDNA and name P"
# Input:
set psf pore6_${sys}_water${water}.psf
set coords run_pore6_${sys}_water${water}.pdb
# Output:
set restFilePrefix smd_${sys}

mol load psf $psf
mol addfile $coords waitfor all
set selAll [atomselect top all]

$selAll set occupancy 0.0
$selAll set beta 0.0

# Select the silicon nitride.
set sel [atomselect top $selText]
$sel set occupancy $occupancy
$sel set beta $occupancy

# Write the smd file.
$selAll writepdb ${restFilePrefix}.pdb

$sel delete
$selAll delete

mol delete top





