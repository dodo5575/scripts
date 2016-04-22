# Cut a double-cone pore in a membrane.
# Use with: tclsh drillPore.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

# Parameters:
set radiusMin 8
set radiusMax 15
# Input:
set pdbIn membrane_hex.pdb
# Output:
set pdbOut pore.pdb

# Get cellBasisVector3_z from the boundary file.
proc readLz {boundaryFile} {
	set in [open $boundaryFile r]
	foreach line [split [read $in] \n] {
    		if {[string match "cellBasisVector3 *" $line]} {
			set lz [lindex $line 3]
			break
		}
	}
	close $in
	return $lz
}

# Determine whether the position {x y z} is inside of pore and
# should be deleted.
proc insidePore {x y z sMin sMax lz} {
	# Get the radius for the double cone at this z-value.
	set s [expr $sMin + 2.0*($sMax-$sMin)/$lz*abs($z)]
	
	return [expr $x*$x + $y*$y < $s*$s]
}

proc drillPore {sMin sMax lz pdbIn pdbOut} {
	set sqrt3 [expr sqrt(3.0)]
	
	# Open the pdb to extract the atom records.
	set out [open $pdbOut w]
	set in [open $pdbIn r]
	set atom 1
	foreach line [split [read $in] \n] {
    		set string0 [string range $line 0 3]
		
		# Just write any line that isn't an atom record.
		if {![string match $string0 "ATOM"]} {
			puts $out $line
			continue
		}
		
		# Extract the relevant pdb fields.
		set serial [string range $line 6 10]
		set x [string range $line 30 37]
		set y [string range $line 38 45]
		set z [string range $line 46 53]
						
		# If atom is outside the pore, write it to the output pdb.
		# Otherwise, exclude it from the resultant pdb.
		if {![insidePore $x $y $z $sMin $sMax $lz]} {
			# Make the atom serial number accurate if necessary.
			if {[string is integer [string trim $serial]]} {
				puts -nonewline $out "ATOM  "
				puts -nonewline $out \
				[string range [format "     %5i " $atom] end-5 end]
				puts $out [string range $line 12 end]
				
			} else {
				puts $out $line
			}
			
			incr atom
		}
	}
	close $in
	close $out
}

set lz [readLz $boundaryFile]
drillPore $radiusMin $radiusMax $lz $pdbIn $pdbOut



