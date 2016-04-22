# Calculate the distance between two CoMs for a trajectory.
# Writes step time(ns) and position (nm) to a file.
# to use: vmd -dispdev text -e stretchDcd.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set dcdFreq 5000
set timestep 1.0; # in femtoseconds
set stride 1
set cutoff 3; # in angstroms
set segA HAIR
set segB HAIR
set resA0 103; # first resid of positive run
set resB0 126; # first resid of negative run
set nBasepairs 8; # last basepair is resA0+nBasepairs-1 and resB0-nBasepairs+1

# Input:
set psf wad_all.psf
set pdb wad_all.pdb
set dcdPrefix /Scr/nanopore/jcomer/myhairpin/race1_dcd/helix_0V
set dcdSuffix ".dcd"
set dcdSet {0}
# Output:
set outPrefix angle_helix_0V0

set degreesPerRadian [expr 180.0/(4.0*atan(1.0))]
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

# Select the atoms defining the bases and backbone.
foreach a $resA b $resB {
    set selBack0($a) [atomselect top "segname $segA and resid $a and name $backAtom"]
    set selBack1($a) [atomselect top "segname $segA and resid [expr $a+1] and name $backAtom"]
    set selBack0($b) [atomselect top "segname $segB and resid $b and name $backAtom"]
    set selBack1($b) [atomselect top "segname $segB and resid [expr $b-1] and name $backAtom"]
    set selBase($a) [atomselect top "segname $segA and resid $a and name $baseAtom"]
    set selBase($b) [atomselect top "segname $segB and resid $b and name $baseAtom"]
}

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
	
	set angle 0.0
	set nAngles 0
	# Determine the angle for bases with at least one H-bond.
	foreach a $resA b $resB {
	    # Count the number of H-bonds.
	    set nBonds 0
	    foreach bA $hb($a) bB $hb($b) {
		set dist [measure bond [list $bA $bB]]
		if {[expr $dist < $cutoff]} {
		    incr nBonds
		}
	    }
	    
	    # Determine the angle.
	    if {$nBonds > 0} {
		incr nAngles
		
		set rBaseA [lindex [$selBase($a) get {x y z}] 0]
		set rBaseB [lindex [$selBase($b) get {x y z}] 0]
		set rBackA0 [lindex [$selBack0($a) get {x y z}] 0]
		set rBackB0 [lindex [$selBack0($b) get {x y z}] 0]
		set rBackA1 [lindex [$selBack1($a) get {x y z}] 0]
		set rBackB1 [lindex [$selBack1($b) get {x y z}] 0]

		set dBack [vecadd [vecsub $rBackA1 $rBackA0] [vecsub $rBackB1 $rBackB0]]
		set dBase [vecsub $rBaseB $rBaseA]
		set comp [expr [vecdot $dBack $dBase]/([veclength $dBack]*[veclength $dBase])]
		set phi [expr $degreesPerRadian*acos($comp)]
		# Keep the angle on the interval [0 90.0]
		if {[expr $phi > 90.0]} {set phi [expr 180.0-$phi]}
		set angle [expr $angle + $phi]
	    }
	}
	
	if {$nAngles > 0} {
	    set angle [expr $angle/$nAngles]

	    # Write the time and angle.
	    puts $out "$t $angle"
	}
	    
	puts -nonewline [format "FRAME %i: " $f]
	puts "$t $angle"
    }
    set nFrames0 [expr $nFrames+$nFrames0]
}

close $out
exit



