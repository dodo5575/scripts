# Add harmonic constraints to silicon nitride.
# to use: vmd -dispdev text -e constrainSilicon.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

# Parameters:
set stride 10
# Spring constant in kcal/(mol A^2)
set occupancy 1.0
set selText "resname DMMP"
# Input:
set psf DMMP_sol.psf
set dcd output/chemx_dmmp_eq2.dcd
# Output:
set outPrefix dmmp

# Load the system.
mol load psf $psf
mol addfile $dcd waitfor all
set nFrames [molinfo top get numframes]
puts [format "Reading %i frames." $nFrames]

# Make the selections.
set selAll [atomselect top all]
set sel [atomselect top $selText]
puts "Performing SMD on [$sel num] atoms."

# Write the smd files.    
for {set f 0} {$f < $nFrames} {incr f $stride} {
    molinfo top set frame $f
    
    $selAll set occupancy 0.0
    $sel set occupancy $occupancy

    # Write the smd file.
    $selAll writepdb ${outPrefix}_f${f}.pdb
}


$sel delete
$selAll delete

mol delete top
exit




