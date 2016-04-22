# Calculate the center of mass for a trajectory.
# to use: vmd -dispdev text -e trackPositionDcd.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set dcdFreq 5000
set selText "resname POT"
set selText1 "resname CLA"
set timestep 1.0
set stride 1
set binSize 5.0

# Input:
set pdb  ../null/02_run/pore+dna-all.pdb
set psf  ../null/02_run/pore+dna-all.psf
set dcdPrefix  /Scr/nanopore/jcomer/myhairpin/null_dcd/run
set dcdSuffix "_6V.dcd"
set dcdSet {10 11}
set xsc ../null/02_run/run3_4V.restart.xsc
# Output:
set outFile ion_null_6V_K.txt
set outFile1 ion_null_6V_Cl.txt

# Read the system size from the xsc file.
# Note: This only works for lattice vectors along the axes!
proc systemSize {xsc} {
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
	
	return [list $lx $ly $lz]
}

proc binTrajectory {histogram binZ stride selText} {
	set bins [llength $histogram]
	set nFrames [molinfo top get numframes]
	
	for {set f 0} {$f < $nFrames} {incr f $stride} {
		molinfo top set frame $f
		
		# Find the density in each bin.
		set z0 [lindex $binZ 0]
		for {set i 1} {$i < $bins} {incr i} {
			set z1 [lindex $binZ $i]
			
			set d0 [lindex $histogram $i]
			set d [getDensity $z0 $z1 $selText]
			#puts "FRAME $f: $d"
			
			lset histogram $i [expr $d+$d0]
		}
		
		puts "FRAME $f:"
		puts $histogram
	}
	
	return $histogram
}

proc getDensity {z0 z1 selText} {
	set sel [atomselect top "($selText) and z >= $z0 and z < $z1"]
	set wat [atomselect top "water and z >= $z0 and z < $z1"]
	return [expr floor(3.0*[$sel num])/[$wat num]]
}

proc main {} {
	global selText selText1 startFrame timestep stride binSize
	global pdb psf dcdPrefix dcdSuffix dcdSet xsc outFile outFile1

	set dcdEnd [llength $dcdSet]	

	# Initialize the histogram bins.
	set l [systemSize $xsc]
	set lz [lindex $l 2]
	set nBins [expr int($lz/$binSize)]
	set nBins2 [expr int($nBins/2)]
	set binZ {}
	set histogram {}
	set histogram1 {}
	for {set i 0} {$i < $nBins} {incr i} {
		lappend binZ [expr ($i-$nBins2-0.5)*$binSize]
		lappend histogram 0
		lappend histogram1 0
	}
	puts "Bin partitions: $binZ\n"

	# Load the system.
	mol load psf $psf pdb $pdb
	set sel [atomselect top $selText]

	# Loop over the dcd files.
	set nFrames0 0
	for {set dcd 0} {$dcd < $dcdEnd} {incr dcd} {
		# Load the trajectory.
		animate delete all
		set dcdNum [lindex $dcdSet $dcd]
		mol addfile "${dcdPrefix}${dcdNum}${dcdSuffix}" waitfor all
		set nFrames [molinfo top get numframes]
		puts [format "Reading %i frames." $nFrames]
		
		set histogram [binTrajectory $histogram $binZ $stride $selText]
		set histogram1 [binTrajectory $histogram1 $binZ $stride $selText1]
		set nFrames0 [expr $nFrames+$nFrames0]
	}
		
	# Normalize the density.
	set h {}
	foreach b $histogram {
		lappend h [expr $b/$nFrames0]
	}
	set h [lrange $h 1 end]
	
	set h1 {}
	foreach b $histogram1 {
		lappend h1 [expr $b/$nFrames0]
	}
	set h1 [lrange $h1 1 end]
		
	set binZ [lrange $binZ 1 end]
	
	# Write the zeroth histogram.
	set out [open $outFile w]
	foreach z $binZ b $h {
		puts -nonewline $out [expr $z-0.5*$binSize]
		puts $out " $b"
	}
	close $out
		
	# Write the first histogram.
	set out [open $outFile1 w]
	foreach z $binZ b $h1 {
		puts -nonewline $out [expr $z-0.5*$binSize]
		puts $out " $b"
	}
	close $out
}

main
exit




