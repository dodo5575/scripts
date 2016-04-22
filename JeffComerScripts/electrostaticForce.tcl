# Calculate the electric force from a grid for a trajectory.
# Writes step time(ns) and force (pN) to a file.
# to use: vmd -dispdev text -e electrostaticForce.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set samplesZ 400
set zStart -200
set zEnd 200
# Input:
set grid open_4V17-20_ghost_volts_ext.dx
switch 2 {
    0 {
	set selText "segname HAIR"
	set psf hairpin_E4V.psf
	set coor coil_first_E4V.pdb
	set outputFile fz_coil_first.dat
    }
    1 {
	set selText "segname HAIR"
	set psf hairpin_E4V.psf
	set coor loop_first_E4V.pdb
	set outputFile fz_loop_first.dat
    }
    2 {
	set selText "segname ADNA BDNA"
	set psf double.psf
	set coor double_E4V.pdb
	set outputFile fz_double.dat
    }
}

source gridForce.tcl
# Read the dx file.
readDx g $grid

# Load the system.
mol load psf $psf
# Load the trajectory.
animate delete all
mol addfile $coor type pdb waitfor all
set sel [atomselect top $selText]
set out [open $outputFile w]
foreach null {0} {set charge [$sel get charge]}
foreach null {0} {set pos [$sel get z]}

# Get the force in piconewtons for this frame.
for {set i 0} {$i < $samplesZ} {incr i} {
    set rz [expr $zStart + ($zEnd-$zStart)/$samplesZ*$i]
    set fz 0.0
    
    foreach r $pos q $charge {
	set rCurr [vecadd $r $rz]
	set ez [interpolateForce g $rCurr 2]
	set fz [expr $fz + $q*$ez]
    }

    set fz [expr 1602.176487*$fz]
    # Write the force in pN and the position in nm.
    puts "electrostatic force: [expr 0.1*$rz] $fz"
    puts $out "[expr 0.1*$rz] $fz"
}

close $out
$sel delete
exit



