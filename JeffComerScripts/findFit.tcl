# This script finds the first frame that fits a molecule in pore.
# to use: vmd -dispdev text -e findFit.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

# Check for atoms of DNA within "distance" of pore atoms.
set distance 2.0
set startFrame 0
set stride 10
set violatorText [format "segname ADNA BDNA and within %.4f of segname SIN" $distance]

# Input:
set porepdb pore1.6.pdb
set porepsf pore1.6.psf
set dnapdb  dsDNA-all.pdb
set dnapsf  dsDNA-all.psf
set dnadcd  slow_ds5.dcd
# Output:
set finalpdb pore1.6+dna.pdb
set finalpsf pore1.6+dna.psf

set temppdb tmp.pdb

# Load the DNA.
mol load psf $dnapsf pdb $dnapdb
set sel [atomselect top "all"]

# Load the DNA trajectory.
animate delete all
mol addfile $dnadcd waitfor all
set nFrames [molinfo top get numframes]
puts [format "Reading %i frames." $nFrames]
set selDNA [atomselect top "segname ADNA BDNA"]
set massCenter [measure center $selDNA weight mass]
puts "mass center of DNA: $massCenter"

# Load the system topology.
package require psfgen
resetpsf
readpsf $porepsf
coordpdb $porepdb
readpsf $dnapsf

# Start at "startFrame" and move forward "stride" steps, looking for the
# first frame that satisfies the constraints.
set n 1
for {set f $startFrame} {$f < $nFrames && $n > 0} {incr f $stride} {
	molinfo top set frame $f
	puts [format "\nFRAME %i" [molinfo top get frame]]
	
	# Load coordinates for this frame.
	$sel writepdb $temppdb
	coordpdb $temppdb	
	writepdb $finalpdb
	writepsf $finalpsf
	
	# Load the system and check for DNA within "distance" of the pore.
	mol load psf $finalpsf pdb $finalpdb
	set violators [atomselect top $violatorText]
	set n [$violators num]
		
	# Delete the combined system.
	mol delete top
	
	puts [format "Frame %i: %i violating atoms" [molinfo top get frame] $n]
}

if {$n != 0} {
	puts "FAILURE: No frame satisfies the constraints!"
	exit
}

# Annouce the selected frame.
puts [format "SUCCESS: Frame %i has been selected!" [molinfo top get frame]]
exit




