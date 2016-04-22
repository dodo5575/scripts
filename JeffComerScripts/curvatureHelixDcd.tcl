# Calculate the distance between two CoMs for a trajectory.
# Writes step time(ns) and position (nm) to a file.
# to use: vmd -dispdev text -e stretchDcd.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set dcdFreq 5000
set selText "all"
set timestep 1.0; # in femtoseconds
set stride 1
set cutoff 3.5; # in angstroms
set segA HAIR
set segB HAIR
set resA0 102; # first resid of positive run
set resB0 127; # first resid of negative run
set nBasepairs 10; # last basepair is resA0+nBasepairs-1 and resB0-nBasepairs+1
set minSamples 5

# Input:
set psf pore+dna-all.psf
set pdb pore+dna-all.pdb
set dcdPrefix  /Projects/alek/HAIRPIN/20_HAIR-PORE/5_RUNS/hairpin-E6V-
set dcdSuffix ".dcd"
set dcdSet {1 2 3 4 5}
# Output:
set outPrefix curveAB_coil_first_6V[lindex $dcdSet 0]-[lindex $dcdSet end]

# Define the hydrogen-bond-forming atoms.
set hydrogenBond {{ADE N1 THY H3} {ADE H61 THY O4} {GUA H21 CYT O2} {GUA H1 CYT N3} {GUA O6 CYT H41}}
# Define the atoms from which the geometry is determined.
set backAtom "C5'"
set baseAtom "C1'"

# Make the hydrogen bond list symmetric by adding reverses.
foreach b $hydrogenBond {
    lappend hydrogenBond [list [lindex $b 2] [lindex $b 3] [lindex $b 0] [lindex $b 1]]
}
# Make list of the first resname.
set hBondResName {}
foreach b $hydrogenBond {
    lappend hBondResName [lindex $b 0]
}

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
set outSamples [open ${outPrefix}.samp w]

# Load the system.
mol load psf $psf pdb $pdb

# Get the residue names.
set resNameA {}
set resNameB {}
foreach a $resA b $resB {
    set selA [atomselect top "segname $segA and resid $a"]
    set selB [atomselect top "segname $segB and resid $b"]
    lappend resNameA [lindex [$selA get resname] 0]
    lappend resNameB [lindex [$selB get resname] 0]
    $selA delete
    $selB delete
}

# Determine the indices of the atoms forming the base hydrogen bonds.
foreach a $resA b $resB rnA $resNameA rnB $resNameB {
    set indList [lsearch -all $hBondResName $rnA]
    
    set indA {}
    set indB {}
    # Find the H-bond forming atoms for these two residues.
    set nBonds 0
    foreach ind $indList {
	foreach {baseA nameA baseB nameB} [lindex $hydrogenBond $ind] {break}
	
	if {[string equal $baseB $rnB]} {
	    set selA [atomselect top "segname $segA and resid $a and name $nameA"]
	    set selB [atomselect top "segname $segB and resid $b and name $nameB"]
	    
	    if {[$selA num] != 1 && [$selB num] != 1}  {
		puts stderr "Warning! Cannot bond $baseA to $baseB."
		puts stderr "Atoms $a:$nameA and $b:$nameB were not found."
	    } else {
		set iA [lindex [$selA get index] 0]
		set iB [lindex [$selB get index] 0]
		lappend indA $iA
		lappend indB $iB
		incr nBonds
		puts "bond: $a $rnA $b $rnB"
	    }
	    $selA delete
	    $selB delete
	}
    }
    set hb($a) $indA
    set hb($b) $indB

    if {$nBonds == 0} {
	puts stderr "Warning! Found no bonds between residues $a $rnA and $b $rnB."
    }
}

# Select the atoms defining the helix shape.
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
	
	# Get the number of H-bonds and positions.
	set nBonds {}
	set posA {}
	set posB {}
	foreach a $resA b $resB {
	    # Count the number of H-bonds on each residue.
	    set n 0
	    foreach bA $hb($a) bB $hb($b) {
		set dist [measure bond [list $bA $bB]]
		if {[expr $dist < $cutoff]} {
		    incr n
		}
	    }
	    lappend nBonds $n
	    lappend posA [measure center $sel($a) weight mass]
	    lappend posB [measure center $sel($b) weight mass]
	}

	# Compute the curvature over the molecule.
	set last [expr $nBasepairs-1]
	set curvature 0.0
	set nUnbroken 0
	for {set j 1} {$j < $last} {incr j} {
	    set jl [expr $j-1]
	    set jr [expr $j+1]
	    # Compute the curvatures for stretches with 3 unbroken bonds.
	    if {[lindex $nBonds $jl]>0 && [lindex $nBonds $jr]>0 && [lindex $nBonds $j]>0} {
		set kA [computeCurvature [lindex $posA $jl] [lindex $posA $j] [lindex $posA $jr]]
		set kB [computeCurvature [lindex $posB $jl] [lindex $posB $j] [lindex $posB $jr]]
		set curvature [expr $curvature + $kA + $kB]
		incr nUnbroken
	    }
	}
	
	if {$nUnbroken > $minSamples} {
	    # Convert to nanometers and compensate for having two strands.
	    set curvature [expr 5.0*$curvature/$nUnbroken]
	    # Write the time (ns) and curvature (nm^-1).
	    puts $out "$t $curvature"
	}
	
	puts $outSamples "$t $nUnbroken"
	puts -nonewline [format "FRAME %i, $nUnbroken samples: " $f]
	puts "$t $curvature"
    }
    set nFrames0 [expr $nFrames+$nFrames0]
}

close $out
close $outSamples
exit



