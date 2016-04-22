# This script will remove all residues in the selection
# from psf and pdf files.
# Use with: vmd -dispdev text -e removeResidues.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set dcdFreq 1000
#set selText "x^2+y^2 < 10^2"
set selText "water"
set timestep 1.0
set stride 1

# Input:
set psf pore1_6_all.psf
set pdb pore1_6_all.pdb
set dcdPrefix /Scr/nanopore/jcomer/myhairpin/pore1.6_open_dcd/heat_p1.6_4V
set dcdSuffix .dcd
set velDcdSuffix .veldcd
set dcdSet {0 1 2 3a 3b 3c 3d 3e 3f 3g}
# Output:
set outFileSuffix open_p1.6_4V[lindex $dcdSet 0]-[lindex $dcdSet end].dat

set dz 20
set z0 -100
set dcdEnd [llength $dcdSet]

set outBin [open bin_water_${outFileSuffix} w]
set out [open centemp_water_${outFileSuffix} w]

# Get the time change between frames in nanoseconds.
set dt [expr 1.0e-6*$timestep*$dcdFreq*$stride]

# Set constants.
set pdbVelFactor 20.45482706
set pdbVelFactor2 [expr $pdbVelFactor*$pdbVelFactor]
set kB 1.3806504e-23
set amu 1.660538782e-27

mol load psf $psf pdb $pdb

# Set up the segments.
set halfDz [expr 0.5*$dz]
set segZ {}
set nz [expr int(ceil(abs(2*$z0)/$dz))+1]
for {set j 0} {$j < $nz} {incr j} {
    lappend segZ [expr $z0 + $j*$dz]
    puts $outBin [expr $z0 + $j*$dz]
}
close $outBin

# Loop over the dcd files.
set nFrames0 0
for {set dcd 0} {$dcd < $dcdEnd} {incr dcd} {

    # Load the position trajectory.
    animate delete all
    set dcdNum [lindex $dcdSet $dcd]
    mol addfile "${dcdPrefix}${dcdNum}${dcdSuffix}" type dcd step $stride waitfor all
    set nFrames [molinfo top get numframes]
    puts [format "Reading %i frames." $nFrames]

    # Load the velocity trajectory.
    mol addfile "${dcdPrefix}${dcdNum}${velDcdSuffix}" type dcd step $stride waitfor all

    # Move forward, computing
    # current at each step.
    for {set f 1} {$f < $nFrames} {incr f} {
	# Get the time in nanoseconds for this frame.
	set t [expr ($nFrames0+$f)*$dt]
	
	# Write the time.
	puts -nonewline $out "$t"

	# Get the temperature in each segment.
	foreach z $segZ {
	    # Make the selection.
	    molinfo top set frame $f
	    set sel [atomselect top "($selText) and abs(z-$z)<$halfDz"]

	    # Get the velocity and other things.
	    molinfo top set frame [expr $f+$nFrames]
	    set vel [$sel get {x y z}]
	    set nAtoms [$sel num]
	    set mass [$sel get mass]
	    $sel delete
	    
	    set eKin 0.0
	    foreach v $vel m $mass {
		set eKin [expr $eKin + 0.5*$m*[veclength2 $v]]
	    }
	    
	    # Convert eKin to joules.
	    set eKin [expr $eKin*$pdbVelFactor2]; # now amu*A^2/ps^2
	    set eKin [expr $eKin*1e-20/1e-24*$amu]
	    
	    # Compute the temperature.
	    set temp [expr 2.0/(3.0*$nAtoms*$kB)*$eKin]
	    puts -nonewline $out " $temp"
	}
	puts $out ""
	puts [format "FRAME %i: $temp" $f]
    }
    set nFrames0 [expr $nFrames+$nFrames0]
}

mol delete top
close $out
exit



