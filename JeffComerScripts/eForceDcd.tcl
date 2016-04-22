# Calculate the electric force from a grid for a trajectory.
# Writes step time(ns) and force (pN) to a file.
# to use: vmd -dispdev text -e stretchDcd.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set dcdFreq 5000
set selText "segname HAIR"
set selTextPos "segname HAIR and resid 102 to 111 118 to 127"
set timestep 1.0
set stride 1

# Input:
set psf nw_pore+dna-all.psf
set pdb nw_pore+dna-all.pdb
set grid /Projects/alek/jcomer/myhairpin/paper/eforce/open_4V17-20_ghost_volts_ext.dx
set dcdPrefix /Scr/nanopore/jcomer/myhairpin/loop_first_dcd/nw_hairpin_first_1.5V
set dcdSuffix ".dcd"
set dcdSet {0 1 2 3 4 5 6 7 8 9}
# Output:
set outPrefix hairpin_first_1.5V[lindex $dcdSet 0]-[lindex $dcdSet end]

source /Projects/alek/jcomer/myhairpin/paper/eforce/vector.tcl
source /Projects/alek/jcomer/myhairpin/paper/eforce/gridForce.tcl
# Read the dx file.
readDx g $grid

# Get the time change between frames in nanoseconds.
set dt [expr 1.0e-6*$timestep*$dcdFreq*$stride]
set dcdEnd [llength $dcdSet]

# Open the output file.
set out [open fz_${outPrefix}.dat w]
set outPosZ [open fzposz_${outPrefix}.dat w]
set outPosS [open fzposs_${outPrefix}.dat w]

# Load the system.
mol load psf $psf pdb $pdb
set sel [atomselect top $selText]
set selPos [atomselect top $selTextPos]
foreach null {0} {set charge [$sel get charge]}

# Loop over the dcd files.
set nFrames0 1
for {set dcd 0} {$dcd < $dcdEnd} {incr dcd} {
    # Load the trajectory.
    animate delete all
    set dcdNum [lindex $dcdSet $dcd]
    mol addfile "${dcdPrefix}${dcdNum}${dcdSuffix}" type dcd step $stride waitfor all
    set nFrames [molinfo top get numframes]
    puts [format "Reading %i frames." $nFrames]

    # Move forward computing the center-of-mass at every step.
    for {set f 0} {$f < $nFrames} {incr f} {
	molinfo top set frame $f

	# Get the time in nanoseconds for this frame.
	set t [expr ($nFrames0+$f)*$dt]
	
	# Get the COM position in nanometers for this frame.
	set r [measure center $selPos weight mass]
	set z [expr 0.1*[lindex $r 2]]
	set s [expr 0.1*sqrt([lindex $r 0]*[lindex $r 0] +[lindex $r 1]*[lindex $r 1])]

	# Get the force in piconewtons for this frame.
	set fz 0
	set pos [$sel get {x y z}]
	foreach r $pos q $charge {
	    set ez [interpolateForce g $r 2]
	    set fz [expr $fz + $q*$ez]
	}
	set fz [expr 1602.176487*$fz]

	# Write the time (ns) and force (pN).
	puts $out "$t $fz"
	# Write the position (nm) and force (pN).
        puts $outPosZ "$z $fz"
	puts $outPosS "$s $fz"

	puts -nonewline [format "FRAME %i: " $f]
	puts "$t $z $fz"
    }
    set nFrames0 [expr $nFrames+$nFrames0]
}
$sel delete
$selPos delete

close $out
close $outPosZ
close $outPosS
exit



