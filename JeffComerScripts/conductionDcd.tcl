# Calculate the current for a trajectory.
# to use: vmd -dispdev text -e conductionDcd.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

# Note: This script is not for general use! It works only for specific systems.
set dcdFreq 5000
set selText "ions"
set selTextDna "segname HAIR"
set timestep 1.0
set nSegments 24
set binZ0 -120.0
set binZ1 120.0

# Input:
set pdb  pore+dna-all.pdb
set psf  pore+dna-all.psf
set dcdPrefix  run
set dcdSuffix "_4V.dcd"
set dcdSet {2 3 4 5 7 8 9}
set xsc run1_4V.restart.xsc
# Output:
set outFile cond_conform_4V.txt
set outFileSeg seg_conform_4V.txt

set dcdEnd [llength $dcdSet]
# Get the time change between frames in femtoseconds.
set dt [expr $timestep*$dcdFreq]

# Prepare the segments.
set binDz [expr ($binZ1-$binZ0)/$nSegments]
set binZ {}
for {set i 0} {$i < $nSegments} {incr i} {
	lappend binZ [expr $binZ0+$i*$binDz] [expr $binZ0+($i+1)*$binDz]
}

proc computeQuantities {} {
	global selDna
	set fourThirdsPi 4.1887902047863
		
	# Compute geometric quantities for the DNA.
	set rc [measure center $selDna]
	foreach {xc yc zc} $rc {break}
	set rgyr [measure rgyr $selDna]
	
	set selDnaIn [atomselect top "segname HAIR and abs(z) < 100"]
	set nDnaIn [$selDnaIn num]
	set rad [$selDnaIn get radius]
	set volDnaIn 0.0
	if {[expr $nDnaIn > 0]} {
		foreach r $rad {
			set volDnaIn [expr $volDnaIn + $fourThirdsPi*($r*$r*$r)]
		}
	}

	# Compute potassium ion quantities.
	set selK [atomselect top "name POT and abs(z) < 100.0"]
	set nK [$selK num]
	set selK [atomselect top "name POT and within 4.0 of segname HAIR"]
	set nK4 [$selK num]
	set selK [atomselect top "name POT and within 8.0 of segname HAIR"]
	set nK8 [$selK num]

	# Compute chlorine ion quantities.
	set selCl [atomselect top "name CLA and abs(z) < 100.0"]
	set nCl [$selCl num]
	set selCl [atomselect top "name CLA and within 4.0 of segname HAIR"]
	set nCl4 [$selCl num]
	set selCl [atomselect top "name CLA and within 8.0 of segname HAIR"]
	set nCl8 [$selCl num]
		
	return [list $xc $yc $zc $rgyr $nDnaIn $volDnaIn $nK $nK4 $nK8 $nCl $nCl4 $nCl8]
}

proc computeSegQuantities {binZ} {
	set resist 0.0
	set quant {}
	# Loop over the segments, computing values for each.
	foreach {z0 z1} $binZ {
		set selDna [atomselect top "segname HAIR and z >= $z0 and z < $z1"]
		set selK [atomselect top "name POT and z >= $z0 and z < $z1"]
		set selCl [atomselect top "name CLA and z >= $z0 and z < $z1"]
		set selH2O [atomselect top "water and z >= $z0 and z < $z1"]
		
		set nDna [$selDna num]
		set nK [$selK num]
		set nCl [$selCl num]
		set nH2O [$selH2O num]
		if {[expr $nK + $nCl != 0]} {
			set resist [expr $resist + 1.0/($nK+$nCl)]
		}
		
		lappend quant $nDna $nK $nCl $nH2O
	}
	
	return [concat $resist $quant]
}

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
set outSeg [open $outFileSeg w]
puts $out "# For $psf with trajectory ${dcdPrefix}(${dcdSet})${dcdSuffix}"
puts $out "# t(ns) I(A) xc yc zc rgyr nDnaIn volDnaIn nK nK4 nK8 nCl nCl4 nCl8"

puts $outSeg "# For $psf with trajectory ${dcdPrefix}(${dcdSet})${dcdSuffix}"
puts $outSeg "# $nSegments segments between $binZ0 and $binZ1"
puts $outSeg "# t(ns) I(A) resist nDna0 nK0 nCl0 nH2O0 nDna1 ..."

# Load the system.
mol load psf $psf pdb $pdb
set sel [atomselect top $selText]
set selDna [atomselect top $selTextDna]

# Loop over the dcd files.
set nFrames0 0
for {set dcd 0} {$dcd < $dcdEnd} {incr dcd} {

# Load the trajectory.
animate delete all
set dcdNum [lindex $dcdSet $dcd]
mol addfile "${dcdPrefix}${dcdNum}${dcdSuffix}" waitfor all
set nFrames [molinfo top get numframes]
puts [format "Reading %i frames." $nFrames]

for {set i 0} {$i < 1} {incr i} {
	# Get the charge of each atom.
	set q [$sel get charge]

	# Get the position data for the first frame.
	molinfo top set frame 0
	set z0 [$sel get z]
}

# Start at frame 1 and move forward, computing
# current at each step.
set n 1
for {set f 1} {$f < $nFrames && $n > 0} {incr f} {
	molinfo top set frame $f
	
	# Get the position data for the current frame.
	set z1 [$sel get z]
	
	# Find the displacements in the z-direction.
	set dz {}
	foreach r0 $z0 r1 $z1 {
		# Compensate for jumps across the periodic cell.
		set z [expr $r1-$r0]
		if {[expr $z > 0.5*$lz]} {set z [expr $z-$lz]}
		if {[expr $z <-0.5*$lz]} {set z [expr $z+$lz]}
		
		lappend dz $z
	}
	
	# Compute the average charge*velocity between the two frames.
	set qvsum [expr [vecdot $dz $q] / $dt]
		
	# We first scale by the system size to obtain the z-current in e/fs.
	set currentZ [expr $qvsum/$lz]
	# Now we convert to amperes.
	set currentZ [expr $currentZ*1.60217733e-4]
	# Get the time in nanoseconds for this frame.
	set t [expr ($nFrames0+$f+1)*$dt*1.e-6]
				
	# Write the current.
	puts "FRAME $f: $t ns"
	puts $out [concat $t $currentZ [computeQuantities]]
	puts $outSeg [concat $t $currentZ [computeSegQuantities $binZ]]
	
	
	# Store the postion data for the next computation.
	set z0 $z1
}
set nFrames0 [expr $nFrames+$nFrames0]
}

close $out
close $outSeg
exit




