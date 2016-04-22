# Calculate the number of charge carriers and current (nA) with multiple
# selections for a trajectory.
# Output format is t(ns) nSel0 ISel0(nA) nSel1 ISel1(nA) ... nTotal ITotal(nA).
# to use: vmd -dispdev text -e multicurrentDcd.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set dcdFreq 5000
set thresh 0.1
set wallText "interpvol0 > ${thresh}"
set backboneAtoms "C1' H1' C2' H2' H2'' C3' O3' H3' C4' O4' H4' C5' O5' H5' H5'' O1P O2P P"
set backboneHydrogen "H1' H2' H2'' H3' H4' H5' H5''"
set selName [list "backbone" "other"]
set selText [list "segname HAIR and name $backboneAtoms" "segname HAIR and not name $backboneAtoms"]
set timestep 1.0

# Input:
set pdb pore+dna_E4V.pdb
set psf pore+dna_E4V.psf
set dxForce pore2.0_nostick1.0_fm.dx
set dcdPrefix /Projects/alek/HAIRPIN/20_HAIR-PORE/5_RUNS/hairpin-
set dcdSuffix ".dcd"
#set dcdSet {E6V-1 E6V-2 E6V-3 E6V-4 E6V-5 E4V-1 E4V-2 E4V-3 E4V-4 E4V-5 E4V-6 E4V-7 E4V-8 E4V-9}
set dcdSet {E6V-3 E6V-4 E6V-5}
# Output:
set outFilePrefix gridforce_${thresh}
set numFilePrefix numforce_${thresh}

proc forceAtoms {selText wallText} {
    set sel [atomselect top "($wallText) and ($selText)"]
    set atoms [$sel get index]
    $sel delete

    return $atoms
}

# Open the output files.
set out {}
set outNum {}
foreach s $selName {
    lappend out [open ${outFilePrefix}_${s}.txt w]
    lappend outNum [open ${numFilePrefix}_${s}.txt w]
}

set dcdEnd [llength $dcdSet]
# Get the time change between frames in femtoseconds.
set dt [expr $timestep*$dcdFreq]

# Load the system.
mol load psf $psf pdb $pdb
mol addFile $dxForce waitfor all

# Loop over the dcd files.
set nFrames0 1
for {set dcd 0} {$dcd < $dcdEnd} {incr dcd} {

    # Load the trajectory.
    animate delete all
    set dcdNum [lindex $dcdSet $dcd]
    mol addfile ${dcdPrefix}${dcdNum}${dcdSuffix} waitfor all
    set nFrames [molinfo top get numframes]
    puts [format "Reading %i frames." $nFrames]

    # Move forward, computing
    # current at each step.
    for {set f 1} {$f < $nFrames} {incr f} {
	molinfo top set frame $f

	# Get the time in nanoseconds for this frame.
	#set t [expr ($nFrames0+$f)*$dt*1.e-6]
	set step [expr $nFrames0+$f]

	# Compute the number of atoms with the applied force.
	foreach s $selText o $out on $outNum {
	    set atoms [forceAtoms $s $wallText]
	    foreach ind $atoms {
		puts $o "$step $ind"
	    }
	    puts $on "$step [llength $atoms]"
	}

	# Write the time.
	puts "FRAME $f: $step [llength $atoms]"
    }
    set nFrames0 [expr $nFrames+$nFrames0]
}

foreach o $out on $outNum {
    close $o
    close $on
}
exit




