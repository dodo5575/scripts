# Bin atom motions in an axially symmetric system.
# to use: vmd -dispdev text -e velocityMapDcd.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set dcdFreq 1000
set selText "name POT"
set timestep 1.0
# The number of radial and axial bins respectively.
set zBins 20
set sBins 20

# Input:
set pdb  poreEmpty-all.pdb
set psf  poreEmpty-all.psf
set dcdPrefix /Scr/alek2/20_DS4/E4V-DS4-noDNA-
set dcdSuffix ".dcd"
set dcdSet {1 2 3 4}
set coords /Scr/alek2/20_DS4/E4V-DS4-noDNA-4.restart.coor
set xsc E4V-DS4-noDNA-4.restart.xsc
# Output:
set outFilePrefix vel_DS4

# Return the maximum radius and z bounds.
proc getExtent {} {
	set all [atomselect top all]
	set minmax [measure minmax $all]
	set r0 [lindex $minmax 0]
	set r1 [lindex $minmax 1]
	foreach {xMin yMin zMin} $r0 {break}
	foreach {xMax yMax zMax} $r1 {break}
	
	if {[expr abs($xMin) > abs($xMax)]} {
		set sx $xMin
	} else {
		set sx $xMax
	}
	if {[expr abs($yMin) > abs($yMax)]} {
		set sy $yMin
	} else {
		set sy $yMax
	}
	
	set sMax [expr sqrt($sx*$sx + $sy*$sy)];
	$all delete
	return [list $sMax $zMin $zMax]
}

proc makeAxialBins {sMax sBins zMin zMax zBins} {
	set dz [expr 1.*($zMax-$zMin)/$zBins]
	set ds [expr 1.*$sMax/$sBins]

	set binS0 {}
	set binS1 {}
	for {set i 0} {$i < $sBins} {incr i} {
		lappend binS0 [expr $i*ds]
		lappend binS1 [expr ($i+1)*ds]
	}
	
	set binZ0 {}
	set binZ1 {}
	for {set i 0} {$i < $zBins} {incr i} {
		lappend binZ0 [expr $i*dz + zMin]
		lappend binZ1 [expr ($i+1)*dz + zMin]	
	}
	
	return [list $binS0 $binS1 $binZ0 $binZ1]
}

proc getCellFromXsc {xsc} {
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
	close $in
	return [list $lx $ly $lz]
}

set dcdEnd [llength $dcdSet]
# Get the time change between frames in femtoseconds.
set dt [expr $timestep*$dcdFreq]

# Load the system.
mol load psf $psf pdb $pdb
mol addfile $coords waitfor all
set sel [atomselect top $selText]

# Get the system size.
set lattice [getCellFromXsc $xsc]
foreach {lx ly lz} $lattice {break}
puts "NOTE: The system size is $lx $ly $lz.\n"

# Prepare the bins.
set ex [getExtent]
foreach {sMax zMin zMax} $ex {break}
set binDz [expr 1.*($zMax-$zMin)/$zBins]
set binDs [expr 1.*$sMax/$sBins]
set sBins1 [expr $sBins-1]
set zBins1 [expr $zBins-1]
	
# Write the bin centers to a file.
# The first data pair in the file is sBins zBins.
# Also make a binZero nested list.
set outBin [open ${outFilePrefix}.bin w]
puts $outBin "$sBins $zBins"
set binZero {}
for {set i 0} {$i < $sBins} {incr i} {
	set sub {}
	set s [expr (0.5+$i)*$binDs]
	
	for {set j 0} {$j < $zBins} {incr j} {
		set z [expr (0.5+$j)*$binDz + $zMin]
		puts $outBin "$s $z"
		
		lappend sub 0
	}
	lappend binZero $sub
}
close $outBin

# Open the output files.
set outVs [open ${outFilePrefix}.vels w]
set outVz [open ${outFilePrefix}.velz w]
set outNum [open ${outFilePrefix}.num w]
set outTime [open "${outFilePrefix}.time" w]

# Loop over the dcd files.
set nFrames0 0
for {set dcd 0} {$dcd < $dcdEnd} {incr dcd} {

# Load the trajectory.
animate delete all
set dcdNum [lindex $dcdSet $dcd]
mol addfile "${dcdPrefix}${dcdNum}${dcdSuffix}" waitfor all
set nFrames [molinfo top get numframes]
puts [format "Reading %i frames." $nFrames]

# Get the positions for the first frame.
foreach i {1} {
	molinfo top set frame 0
	set r0 [$sel get {x y z}]
}

# Move forward, putting the atoms in bins.
for {set f 1} {$f < $nFrames} {incr f} {
	molinfo top set frame $f
	
	# Empty the bins.
	set binVs $binZero
	set binVz $binZero
	set binNum $binZero
		
	# Get the time in nanoseconds for this frame.
	set t [expr ($nFrames0+$f)*$dt*1.e-6]

	# Write the time.
	puts "FRAME $f: $t"
	
	# Get the positions of the selection.
	set r [$sel get {x y z}]
	foreach p $r p0 $r0 {
		set s [expr sqrt([lindex $p 0]*[lindex $p 0]+[lindex $p 1]*[lindex $p 1])]
		set z [lindex $p 2]
		set bs [expr int($s/$binDs)]
		set bz [expr int(($z-$zMin)/$binDz)]
		
		# Bound the indices.
		if {$bs < 0} {set bs 0}
		if {$bs > $sBins1} {set bs $sBins1}
		if {$bz < 0} {set bz 0}
		if {$bz > $zBins1} {set bz $zBins1}
		
		set dr [vecsub $p $p0]
		foreach {dx dy dz} $dr {break}
		# Compensate for jumps across the periodic cell.
		if {[expr $dz > 0.5*$lz]} {set dz [expr $dz-$lz]}
		if {[expr $dz <-0.5*$lz]} {set dz [expr $dz+$lz]}
		
		set vs [expr sqrt($dx*$dx + $dy*$dy)/$dt]
		set vz [expr $dz/$dt]
		
		# Add the radial velocity to this bin.
		set x [lindex $binVs $bs $bz]
		lset binVs $bs $bz [expr $x + $vs]
		
		# Add the axial velocity to this bin.
		set x [lindex $binVz $bs $bz]
		lset binVz $bs $bz [expr $x + $vz]
		
		# Increment the number in this bin.
		set x [lindex $binNum $bs $bz]
		lset binNum $bs $bz [expr $x + 1]
	}
	
	# Save r for the next frame.
	set r0 $r
	
	# Write the time.
	puts $outTime $t
		
	# Write vs.
	foreach j $binVs {
		foreach item $j {
			puts $outVs $item
		}
	}
	
	# Write vz.
	foreach j $binVz {
		foreach item $j {
			puts $outVz $item
		}
	}
	
	# Write num.
	foreach j $binNum {
		foreach item $j {
			puts $outNum $item
		}
	}
}
set nFrames0 [expr $nFrames+$nFrames0]
}

close $outVs
close $outVz
close $outNum
close $outTime
exit



