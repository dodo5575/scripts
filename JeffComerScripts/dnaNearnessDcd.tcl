# Calculate the number of charge carriers and current (nA) with multiple
# selections for a trajectory.
# Output format is t(ns) nSel0 ISel0(nA) nSel1 ISel1(nA) ... nTotal ITotal(nA).
# to use: vmd -dispdev text -e multicurrentDcd.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set dcdFreq 5000
set wallText "resname SIN"
set selName [list "coil" "helix"]
set selText [list "segname HAIR and resid >= 128" "segname HAIR and resid < 128 and not resid 116 115 114 113"]
set radius0 {0.0 2.0 2.2 2.4 2.6 2.8}
set radius1 {2.0 2.2 2.4 2.6 2.8 3.0}
set timestep 1.0

# Input:
set pdb pore+dna-all.pdb
set psf pore+dna-all.psf
set dcdPrefix /Projects/alek/HAIRPIN/20_HAIR-PORE/5_RUNS/hairpin-
set dcdSuffix ".dcd"
set dcdSet {E6V-3 E6V-4 E6V-5}
#set dcdSet {E6V-4}
set xsc hairpin-E4V-6.restart.xsc
# Output:
set outFilePrefix stick_4-6V

proc numAtoms {selText wallText r0 r1} {
    set sel [atomselect top "($selText) and (within $r1 of $wallText) and not (within $r0 of $wallText)"]
    set n [$sel num]
    $sel delete

    return $n
}

# Write the radii to a file.
set outRad [open $outFilePrefix.rad w]
puts $outRad $radius0
puts $outRad $radius1
close $outRad

# Open the output files.
set out {}
foreach s $selName {
    lappend out [open ${outFilePrefix}_${s}.txt w]
}

set dcdEnd [llength $dcdSet]
# Get the time change between frames in femtoseconds.
set dt [expr $timestep*$dcdFreq]

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

# Move forward, computing
# current at each step.
for {set f 1} {$f < $nFrames} {incr f} {
    molinfo top set frame $f

    # Get the time in nanoseconds for this frame.
    set t [expr ($nFrames0+$f)*$dt*1.e-6]

    # Write the time.
    puts "FRAME $f: $t"

    # Compute the density of DNA at distances.
    foreach s $selText o $out {
	puts -nonewline $o $t
	
	foreach r0 $radius0 r1 $radius1 {
	    set n [numAtoms $s $wallText $r0 $r1]
	    
	    puts -nonewline $o " $n"
	    puts -nonewline " $n"
	}
	puts $o ""
	puts ""
    }
}
set nFrames0 [expr $nFrames+$nFrames0]
}

foreach o $out {
    close $o
}
exit




