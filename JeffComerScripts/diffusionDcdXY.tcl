# Calculate the current for a trajectory.
# Writes time(ns) current(nA) to a file.
# to use: vmd -dispdev text -e currentDcd.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set dcdFreq 5000
set selText "name OH2"
set timestep 1.0
set stride 10

set nameList {open_p2.0_4V}
set name [lindex $nameList 0]
# Input:
set psf pore2.0_all.psf
set pdb pore2.0_all.pdb
set dcdPrefix /Scr/nanopore/jcomer/myhairpin/pore2.0_open_dcd/$name
set dcdSuffix .dcd
set dcdSet {2}
set xsc output/eq3.restart.xsc
# Output:
set outFile diffxy${stride}_${name}[lindex $dcdSet 0]-[lindex $dcdSet end].dat

set dcdEnd [llength $dcdSet]
# Get the time change between frames in nanoseconds.
set dt [expr 1.0e-6*$timestep*$dcdFreq*$stride]

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

# Open the output file.
set out [open $outFile w]

# Load the system.
mol load psf $psf pdb $pdb
set sel [atomselect top $selText]

# Loop over the dcd files.
set nFrames0 1
for {set dcd 0} {$dcd < $dcdEnd} {incr dcd} {
    # Load the trajectory.
    animate delete all
    set dcdNum [lindex $dcdSet $dcd]
    mol addfile "${dcdPrefix}${dcdNum}${dcdSuffix}" type dcd step $stride waitfor all
    set nFrames [molinfo top get numframes]
    puts [format "Reading %i frames." $nFrames]

    # Get the position data for the first frame.
    molinfo top set frame 0
    foreach zero {0} {set r0 [$sel get {x y z}]}

    # Move forward, computing
    # current at each step.
    for {set f 1} {$f < $nFrames} {incr f} {
	molinfo top set frame $f
	# Get the time in nanoseconds for this frame.
	set t [expr ($nFrames0+$f-0.5)*$dt]

	# Get the position data for the current frame.
	set r1 [$sel get {x y z}]
	
	# Find the displacements in each direction.
	set rSqSum 0.0
	foreach p0 $r0 p1 $r1 {
	    set d [vecsub $p1 $p0]
	    foreach {dx dy dz} $d {break}
	    
	    # Handle the periodic boundaries.
	    if {[expr $dx > 0.5*$lx]} {set dx [expr $dx-$lx]}
	    if {[expr $dx <-0.5*$lx]} {set dx [expr $dx+$lx]}
	    if {[expr $dy > 0.5*$ly]} {set dy [expr $dy-$ly]}
	    if {[expr $dy <-0.5*$ly]} {set dy [expr $dy+$ly]}
	    if {[expr $dz > 0.5*$lz]} {set dz [expr $dz-$lz]}
	    if {[expr $dz <-0.5*$lz]} {set dz [expr $dz+$lz]}
	    
	    set rSqSum [expr $rSqSum + $dx*$dx + $dy*$dy]
	}
	set rSq [expr $rSqSum/[llength $r0]]
	
	set dif [expr 0.01*$rSq/(4.0*$dt)]
	
	# Write the current.
	puts $out "$t $dif"
	puts -nonewline [format "FRAME %i: " $f]
	puts "$t $dif"
	
	# Store the position data for the next computation.
	set r0 $r1
    }
    set nFrames0 [expr $nFrames+$nFrames0]
}

$sel delete
mol delete top
close $out
exit




