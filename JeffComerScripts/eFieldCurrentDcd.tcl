# Compute the electric dipole moment for a trajectory.
# to use: vmd -dispdev text -e dipoleDcd.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set dcdFreq 5000
set selText "name POT and not within 6 of resname SIN"
set muK 6.071e-8
set muCl 6.678e-8
set mu $muK
set e 1.60217653e-19;
set startFrame 0
set timestep 1.0
set stride 1

# Input:
set pdb  poreEmpty-all.pdb
set psf  poreEmpty-all.psf
set dcdPrefix /Scr/nanopore/jcomer/myhairpin/open_dcd/run
set dcdSuffix "_4V.dcd"
set dcdSet {9 10 11 12}
set xsc run12_4V.restart.xsc
set eFieldFile open_10-12_z.eField
# Output:
set outFile ecurrK_open_[lindex $dcdSet 0]-[lindex $dcdSet end].txt

# Get the time change between frames in femtoseconds.
set dt [expr $timestep*$dcdFreq]
set dcdEnd [llength $dcdSet]

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

# Read the z-electric field as a function of z, Ez(z).
set in [open $eFieldFile r]
set binZ {}
set eField {}
foreach line [split [read $in] "\n"] {
	if {[llength $line] > 1} {
		set param [split $line]
		lappend binZ [lindex $param 0]
		lappend eField [lindex $param 1]
	}
}
puts "NOTE: Read [llength $binZ] electric field values."
close $in

# Organize the eField bins.
set z0 [lindex $binZ 0]
set dz [expr [lindex $binZ 1] - $z0]
set nBins [llength $binZ]
puts "Bin width in angstroms: $dz"

# Open the output file.
set out [open $outFile w]

# Load the system.
mol load psf $psf pdb $pdb

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
for {set f $startFrame} {$f < $nFrames} {incr f $stride} {
	molinfo top set frame $f

	set curr 0.0
	set sel [atomselect top $selText]
	set posZ [$sel get z]
	$sel delete
	
	foreach z $posZ {
		# Determine which bin we are in.
		set bin  [expr int(($z-$z0)/$dz + 0.5)]
		if {[expr $bin < 0 || $bin > $nBins-1]} {
			continue
		}
	
		# Add this contribution to the current.
		set curr [expr $curr + [lindex $eField $bin]]
	}
	
	# Multiply constants to get the result in nA.
	set curr [expr $curr*$mu*$e*1.e29/$lz]
		
	# Get the time in nanoseconds for this frame.
	set t [expr ($nFrames0+$f)*$dt*1.e-6]
	
	puts $out "$t $curr"
	puts -nonewline [format "FRAME %i: " $f]
	puts "$t $curr"
	
}
set nFrames0 [expr $nFrames+$nFrames0]
}

close $out
exit




