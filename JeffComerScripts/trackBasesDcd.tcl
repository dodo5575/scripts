# Calculate the center of mass for a trajectory.
# Writes `step time(ns) isPaired angle(rad) pairDist(A)' to a file.
# to use: vmd -dispdev text -e trackBasesDcd.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set dcdFreq 5000
set selText "segname HAIR"
set timestep 1.0
set stride 1
set cutoff 3.0
set resIdList0 {102 103 104 105 106 107 108 109 110 111}
set resIdList1 {127 126 125 124 123 122 121 120 119 118}
set up {0.0 0.0 1.0}

# Input:
set pdb  tall_sys0.pdb
set psf  tall_sys0.psf
set dcdPrefix uz0_zip
set dcdSuffix ".dcd"
set dcdSet {0 1 2 3}
# Output:
set outPrefix base_uz0

set bondA "name N1"
set bondT "name H3"
set bondC "name N3"
set bondG "name H1"

# Get the time change between frames in femtoseconds.
set dt [expr $timestep*$dcdFreq]
set dcdEnd [llength $dcdSet]

# Open the output file.
set out {}
foreach res1 $resIdList1 {
    lappend out [open ${outPrefix}_${res1}.dat w]
}

# Load the system.
mol load psf $psf pdb $pdb

# Get the residue names.
set resNameList0 {}
set resNameList1 {}
foreach res0 $resIdList0 res1 $resIdList1 {
    set sel [atomselect top "($selText) and resid $res0"]
    lappend resNameList0 [lindex [$sel get resname] 0]
    $sel delete

    set sel [atomselect top "($selText) and resid $res1"]
    lappend resNameList1 [lindex [$sel get resname] 0]
    $sel delete
}
puts "resNameList0: $resNameList0"
puts "resNameList1: $resNameList1"

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
	
	# Get the time in nanoseconds for this frame.
	set step [expr ($nFrames0+$f)*$dcdFreq]
	set t [expr ($nFrames0+$f)*$dt*1.e-6]

	set outText {}
	lappend $outText $t
	set angle0 0.0
	set angle1 0.0
	set numPaired 0

	foreach res0 $resIdList0 res1 $resIdList1 base0 $resNameList0 o $out {
	    set text0 "($selText) and resid $res0"
	    set text1 "($selText) and resid $res1"

	    # Determine whether the bases are paired.
	    if {[string equal $base0 ADE]} {
		set pair0 "$text0 and $bondA"
		set pair1 "$text1 and $bondT"
		set connect0 N9
		set connect1 N1
	    } elseif {[string equal $base0 CYT]} {
		set pair0 "$text0 and $bondC"
		set pair1 "$text1 and $bondG"
		set connect0 N1
		set connect1 N9
	    } elseif {[string equal $base0 GUA]} {
		set pair0 "$text0 and $bondG"
		set pair1 "$text1 and $bondC"
		set connect0 N9
		set connect1 N1
	    } elseif {[string equal $base0 THY]} {
		set pair0 "$text0 and $bondT"
		set pair1 "$text1 and $bondA"
		set connect0 N1
		set connect1 N9
	    } else {
		puts "ERROR: Unrecognized base $base0"
		exit
	    }

	    set selA [atomselect top $pair0]
	    set selB [atomselect top $pair1]
	    set rA [lindex [$selA get {x y z}] 0]
	    set rB [lindex [$selB get {x y z}] 0]
	    set pairDist [veclength [vecsub $rB $rA]]
	    set paired [expr $pairDist < $cutoff]
	    $selA delete
	    $selB delete

	    # Determine the angle of the bases.
	    set selA [atomselect top "$text1 and name C1'"]
	    set selB [atomselect top "$text1 and name $connect1"]
	    set rA [lindex [$selA get {x y z}] 0]
	    set rB [lindex [$selB get {x y z}] 0]
	    set rBase1 [vecsub $rB $rA]
	    set comp1 [vecscale [expr 1./[veclength $rBase1]] [vecdot $rBase1 $up]]
	    set angle1 [expr acos($comp1)]
	    $selA delete
	    $selB delete

	    # Write the results.
	    puts $o "$step $t $paired $angle1 $pairDist"

	    set numPaired [expr $numPaired + $paired]
	}
	
	puts -nonewline [format "FRAME %i: " $f]
	puts "$step $numPaired"
    }
    set nFrames0 [expr $nFrames+$nFrames0]
}

foreach o $out {
    close $o
}

exit




