# Compute the concentration of solute in a zone for a trajectory.
# Output time(ns) concentration(mol/kg).
# to use: vmd -dispdev text -e concDcd.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set dcdFreq 5000
set zoneText "z < -135 and z > -155"
set soluteText "name POT CLA"
set waterText "resname TIP3 and name OH2"
set timestep 1.0
set stride 1

# Input:
set pdb  poreEmpty-all.pdb
set psf  poreEmpty-all.psf
set dcdPrefix run
set dcdSuffix "_4V.dcd"
set dcdSet {5 6 7 8 9 10 11 12}
# Output:
set outFile conc_open.txt

# Get the time change between frames in femtoseconds.
set dt [expr $timestep*$dcdFreq]
set dcdEnd [llength $dcdSet]

# Open the output file.
set out [open $outFile w]

# Load the system.
mol load psf $psf pdb $pdb

# Loop over the dcd files.
set nFrames0 1
for {set dcd 0} {$dcd < $dcdEnd} {incr dcd} {

# Load the trajectory.
animate delete all
set dcdNum [lindex $dcdSet $dcd]
mol addfile "${dcdPrefix}${dcdNum}${dcdSuffix}" waitfor all
set nFrames [molinfo top get numframes]
puts [format "Reading %i frames." $nFrames]

# Start at "startFrame" and move forward, computing
# the sasa at each step.
for {set f 0} {$f < $nFrames} {incr f $stride} {
	molinfo top set frame $f

	# Get the number of molecules of solute and water in the zone.
	set sel [atomselect top "(${zoneText}) and (${soluteText})"]
	set water [atomselect top "(${zoneText}) and (${waterText})"]
	set ns [$sel num]
	set nw [$water num]
	$sel delete
	$water delete
	
	# Compute the concentration in mol/kg.
	set a [expr 55.523*$ns/$nw]
		
	# Get the time in nanoseconds for this frame.
	set t [expr ($nFrames0+$f)*$dt*1.e-6]
	
	puts $out "$t $a"
	puts -nonewline [format "FRAME %i: " $f]
	puts "$t $a"
	
}
set nFrames0 [expr $nFrames+$nFrames0]
}

close $out
exit




