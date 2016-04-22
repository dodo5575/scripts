# Calculate the center of mass for a trajectory.
# Writes step time(ns) and position (nm) to a file.
# to use: vmd -dispdev text -e trackPositionDcd.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set type OSI
set bonds 2
set dcdFreq 5000
set selText "(type ${type} and numbonds < $bonds) and within 4.0 of water"
set timestep 1.0
set stride 1

# Input:
set psf dmmp+raw.psf
set pdb dmmp+raw.pdb
set dcdPrefix output/dmmp_raw_eq
set dcdSuffix ".dcd"
set dcdSet {2}
# Output:
set outPrefix dangling_${type}

# Get the time change between frames in nanoseconds.
set dt [expr 1.0e-6*$timestep*$dcdFreq*$stride]
set dcdEnd [llength $dcdSet]

# Open the output file.
set out [open $outPrefix.dat w]

# Load the system.
mol load psf $psf pdb $pdb
set pickList {}

# Loop over the dcd files.
set nFrames0 0
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

	# Get the number in the selection.
	set sel [atomselect top $selText]
	foreach index [$sel get index] {
	    lappend pickList $index
	}
	$sel delete

	set pickList [lsort -integer -unique $pickList]
	set n [llength $pickList]

	# Write the time and position.
	puts $out "$t $n"
	
	puts -nonewline [format "FRAME %i: " $f]
	puts "$t $n"
    }
    set nFrames0 [expr $nFrames+$nFrames0]
}
close $out

# Write the picked atoms.
set outPick [open $outPrefix.pick w]
puts $outPick $pickList
close $outPick

exit
