# Calculate the number of charge carriers and current (nA) with multiple
# selections for a trajectory.
# Output format is t(ns) nSel0 ISel0(nA) nSel1 ISel1(nA) ... nTotal ITotal(nA).
# to use: vmd -dispdev text -e multicurrentDcd.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set dcdFreq 5000
set wallText "resname SIN"
set selName [list "K" "Cl"]
set selText [list "name POT" "name CLA"]
set timestep 1.0
set dz 20.0
set z0 -120.0
set radius0 [expr 12.1 - 3.0]
set radius1 [expr $radius0 + 17.63]
set lengthPore 200.0

# Input:
set pdb poreEmpty-all.pdb
set psf poreEmpty-all.psf
set dcdPrefix /Scr/nanopore/jcomer/myhairpin/open_dcd/run
set dcdSuffix ".dcd"
set dcdSet {5_4V 6_4V 7_4V 8_4V 9_4V 10_4V 11_4V 12_4V 13 14 15 16 17 18 19 20 21}
#set dcdSet {5_4V}
set xsc run21.restart.xsc
# Output:
set outFilePrefix snum_open

# Read the system size from the xsc file.
# Note: This only works for lattice vectors along the axes!
set in [open $xsc r]
foreach line [split [read $in] "\n"] {
	if {![string match "#*" $line]} {
		set param [split $line]
		puts $param
		set lx [lindex $param 1]
		set ly [lindex $param 5]
		set lz [lindex $param 9]
		break
	}
}
puts "NOTE: The system size is $lx $ly $lz.\n"
close $in
set nSegments [expr int(ceil((0.5*$lz-$z0)/$dz))]
set slope [expr 2.0*($radius1-$radius0)/$lengthPore]

# Write the bins.
set outBin [open ${outFilePrefix}_rad.txt w]
for {set i 0} {$i < $nSegments} {incr i} {
    set zi [expr $dz*$i + $z0]
    set zi1 [expr $dz*($i+1) + $z0]
    puts $outBin "$zi $zi1"
}
close $outBin

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
set sel {}
foreach s $selText {
    lappend sel [atomselect top $s]
}

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
	foreach s $sel o $out {
	    # Zero the segments.
	    for {set i 0} {$i < $nSegments} {incr i} {
		set hist($i) 0
	    }
    
	    set pos [$s get {x y z}]
	    foreach r $pos {
		foreach {x y z} $r {break}
		set s [expr sqrt($x*$x + $y*$y)]
		set sPore [expr $radius0 + $slope*abs($z)]
		
		set seg [expr int(floor(($z-$z0)/$dz))]
		if {$s > $sPore && $seg >= 0 && $seg < $nSegments} {
		    incr hist($seg)
		}
	    }
	    
	    # Write the results.
	    puts -nonewline $o "$t"
	    for {set i 0} {$i < $nSegments} {incr i} {
		puts -nonewline $o " $hist($i)"
	    }
	    puts $o ""
	}
    }
    set nFrames0 [expr $nFrames+$nFrames0]
}

foreach o $out {
    close $o
}
foreach s $sel {
    $s delete
}
exit



