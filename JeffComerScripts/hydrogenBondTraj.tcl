# Calculate the distance between two CoMs for a trajectory.
# Writes step time(ns) position (nm) cross-sectional radius (nm) to a file.
# to use: vmd -dispdev text -e stretchDcd.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

# Requires: stride timestep dcdFreq psf pdb xsc dcd outSuffix moiety displayPeriod
set cutoff 2.6
#set cutoff 4

switch $moiety {
    hpDNA {
	set segA HAIR
	set segB HAIR
	set resA0 102; # first resid of positive run
	set resB0 127; # first resid of negative run
	set nBasepairs 10; # last basepair is resA0+nBasepairs-1 and resB0-nBasepairs+1
    }
    dsDNA {
	set segA ADNA
	set segB BDNA
	set resA0 1; # first resid of positive run
	set resB0 111; # first resid of negative run
	set nBasepairs 10; # last basepair is resA0+nBasepairs-1 and resB0-nBasepairs+1	
    }
    default {
	puts stderr "ERROR: Unrecognized moiety $moiety."
	exit
    }
}

# Define the hydrogen-bond-forming atoms.
set hydrogenBond {{ADE N1 THY H3} {ADE H61 THY O4} {GUA H21 CYT O2} {GUA H1 CYT N3} {GUA O6 CYT H41}}

# Open the output file.
set out [open output/hbond_${outSuffix} w]
set outBroken [open output/broken_hbond_${outSuffix} w]
set outIndex [open output/index_hbond_${outSuffix} w]

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
set hbA {}
set hbB {}
foreach a $resA b $resB rnA $resNameA rnB $resNameB {
    set indList [lsearch -all $hBondResName $rnA]
    
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
		lappend hbA $iA
		lappend hbB $iB
		incr nBonds
		puts "bond: $a $rnA $b $rnB"
		puts $outIndex "$iA $a $rnA $iB $b $rnB"
	    }
	    $selA delete
	    $selB delete
	}
    }

    if {$nBonds == 0} {
	puts stderr "Warning! Found no bonds between residues $a $rnA and $b $rnB."
    }
}
puts "Found [llength $hbA] potential hydrogen bonds!\n"
close $outIndex

# Loop over the dcd files.
set nFrames0 0
foreach dcdFile $dcd {
    # Load the trajectory.
    animate delete all
    mol addfile $dcdFile type dcd step $stride waitfor all
    set nFrames [molinfo top get numframes]
    puts [format "Reading %i frames." $nFrames]

    # Move forward computing the center-of-mass at every step.
    for {set f 0} {$f < $nFrames} {incr f} {
	molinfo top set frame $f

	# Get the time in nanoseconds for this frame.
	set t [expr ($nFrames0+$f)*$dt]

	puts -nonewline $outBroken $t
	# Determine the number of H-bonds.
	set n 0
	foreach a $hbA b $hbB {
	    set dist [measure bond [list $a $b]]
	    if {[expr $dist < $cutoff]} {
		incr n
	    } else {
		puts -nonewline $outBroken " $a $b"
	    }
	}
	puts $outBroken ""

	# Write the time and distance.
	puts $out "$t $n"
	if {$f % $displayPeriod == 0} {
	    puts -nonewline [format "FRAME %i: " $f]
	    puts "$t $n"
	}
    }
    set nFrames0 [expr $nFrames+$nFrames0]
}

mol delete top
close $out
close $outBroken



