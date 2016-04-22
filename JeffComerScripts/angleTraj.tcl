# Calculate the distance between two CoMs for a trajectory.
# Writes step time(ns) position (nm) cross-sectional radius (nm) to a file.
# to use: vmd -dispdev text -e stretchDcd.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

# Requires: stride timestep dcdFreq psf pdb dcd outSuffix displayPeriod

#set selText "segname HAIR and resid 102 to 111 118 to 127 and backbone"
#set selText "segname HAIR and resid 102 to 111 118 to 127"
set selText "name N1 C2 N3 C4 C5 C6" 
set segA HAIR
set segB HAIR
set startA 102
set startB 127
set nPairs 10
set psf0 wad_all.psf
set pdb0 helix_0V0.pdb

# Get the time change between frames in nanoseconds.
set dt [expr 1.0e-6*$timestep*$dcdFreq*$stride]
set dcdEnd [llength $dcd]
set pi [expr 4.0*atan(1.0)]

# Open the output file.
set out [open angle/angle_${outSuffix} w]

# Load the system.
mol load psf $psf pdb $pdb
set endA [expr $startA + $nPairs - 1]
set endB [expr $startB - $nPairs + 1]
set sel0 [atomselect top "backbone and ((segname $segA and resid $startA) or (segname $segB and resid $startB))"]
set sel1 [atomselect top "backbone and ((segname $segA and resid $endA) or (segname $segB and resid $endB))"]

puts "Start: [$sel0 num] atoms"
puts "End: [$sel1 num] atoms"

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

	# Get the helix direction.
	set r0 [measure center $sel0 weight mass]
	set r1 [measure center $sel1 weight mass]
	set dir [vecsub $r1 $r0]
	set dir [vecscale [expr 1.0/[veclength $dir]] $dir]

	set angle 0.0
	for {set i 0} {$i < $nPairs} {incr i} {
	    set resA [expr $startA + $i]
	    set resB [expr $startB - $i]
	    set selA [atomselect top "segname $segA and resid $resA and ($selText)"]
	    set selB [atomselect top "segname $segB and resid $resB and ($selText)"]

	    set rA [measure center $selA]
	    set rB [measure center $selB]
	    $selA delete
	    $selB delete

	    set d [vecsub $rB $rA]
	    set dz [vecdot $d $dir]
	    set dm [veclength $d]
	    set phi [expr asin($dz/$dm)]
	    set angle [expr $angle + $phi]
	}

	# Write the angle.
	set angleMean [expr $angle*180.0/$pi/$nPairs]
	puts $out "$t $angleMean"
		
	if {$f % $displayPeriod == 0} {
	    puts -nonewline [format "FRAME %i: " $f]
	    puts "$t $angleMean"
	}
    }
    set nFrames0 [expr $nFrames+$nFrames0]
}

$sel0 delete
$sel1 delete
mol delete top
close $out



