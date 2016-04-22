# Calculate the center of mass for a trajectory.
# Writes time(ns) position(nm) to a file.
# to use: vmd -dispdev text -e trackPositionDcd.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set dcdFreq 5000
set selText "segname HAIR"
set timestep 1.0
set stride 1
set cutoff 3.0

# Input:
set pdb  tall_sys1.pdb
set psf  tall_sys1.psf
set dcdPrefix uz1_zip
set dcdSuffix ".dcd"
set dcdSet {0 1 2 3}
# Output:
set outFile pairs_${dcdPrefix}.txt

set bondA "resname ADE and name N1"
set bondT "resname THY and name H3"
set bondC "resname CYT and name N3"
set bondG "resname GUA and name H1"
set pairAT "$selText and $bondA and within $cutoff of ($bondT)"
set pairCG "$selText and $bondC and within $cutoff of ($bondG)"

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
	
	set numPairs 0
	set selAT [atomselect top $pairAT]
	set selCG [atomselect top $pairCG]
	
	set numPairs [expr $numPairs + [$selAT num]]
	set numPairs [expr $numPairs + [$selCG num]]
	$selAT delete
	$selCG delete
	
	# Get the time in nanoseconds for this frame.
	set t [expr ($nFrames0+$f)*$dt*1.e-6]
	
	puts $out "$t $numPairs"
	puts -nonewline [format "FRAME %i: " $f]
	puts "$t $numPairs"
}
set nFrames0 [expr $nFrames+$nFrames0]
}

close $out
exit




