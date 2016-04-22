# Calculate the distance between two CoMs for a trajectory.
# Writes step time(ns) and position (nm) to a file.
# to use: vmd -dispdev text -e stretchDcd.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set dcdFreq 5000
set selText "all"
set timestep 1.0; # in femtoseconds
set stride 10
set segA HAIR
set segB HAIR
set resA0 102; # first resid of positive run
set resB0 127; # first resid of negative run
set nBasepairs 10; # last basepair is resA0+nBasepairs-1 and resB0-nBasepairs+1

# Input:
set psf pore+dna-all.psf
set pdb pore+dna-all.pdb
set dcdPrefix  /Projects/alek/HAIRPIN/20_HAIR-PORE/5_RUNS/hairpin-E6V-
set dcdSuffix ".dcd"
set dcdSet {1 2 3 4 5}
# Output:
set outPrefix curvall_coil_first_6V[lindex $dcdSet 0]-[lindex $dcdSet end]

# Find the complementary bases.
set resA {}
set resB {}
for {set i 0} {$i < $nBasepairs} {incr i} {
    lappend resA [expr $resA0 + $i]
    lappend resB [expr $resB0 - $i]
}

# Get the time change between frames in nanoseconds.
set dt [expr 1.0e-6*$timestep*$dcdFreq*$stride]
set dcdEnd [llength $dcdSet]

# Open the output files.
set out [open ${outPrefix}.dat w]

# Load the system.
mol load psf $psf pdb $pdb

# Select the atoms defining the molecule.
foreach a $resA b $resB {
    set sel($a) [atomselect top "segname $segA and resid $a and ($selText)"]
    set sel($b) [atomselect top "segname $segB and resid $b and ($selText)"]
}

proc computeCurvature {r0 r1 r2} {
    set rl [vecsub $r0 $r1]
    set rr [vecsub $r2 $r1]
    set d [vecadd $rl $rr]
    set dUnit [vecscale [expr 1.0/[veclength $d]] $d]
    return [expr 2.0*[vecdot $rr $dUnit]/[veclength2 $rr]]
}

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
	
	# Get the positions.
	set posA {}
	set posB {}
	foreach a $resA b $resB {
	    lappend posA [measure center $sel($a) weight mass]
	    lappend posB [measure center $sel($b) weight mass]
	}

	# Compute the curvature over the molecule.
	set n 0
	set last [expr $nBasepairs-1]
	set curvature 0.0
	for {set j 1} {$j < $last} {incr j} {
	    set jl [expr $j-1]
	    set jr [expr $j+1]
	    set kA [computeCurvature [lindex $posA $jl] [lindex $posA $j] [lindex $posA $jr]]
	    set kB [computeCurvature [lindex $posB $jl] [lindex $posB $j] [lindex $posB $jr]]
	    set curvature [expr $curvature + $kA + $kB]
	    incr n
	}
       
	# Convert to nanometers^-1 and compensate for having two strands.
	set curvature [expr 5.0*$curvature/$n]
	# Write the time (ns) and curvature (nm^-1).
	puts $out "$t $curvature"
		
	puts -nonewline [format "FRAME %i: " $f]
	puts "$t $curvature"
    }
    set nFrames0 [expr $nFrames+$nFrames0]
}

close $out
exit



