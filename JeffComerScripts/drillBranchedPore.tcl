# Cut a branched pore in a membrane.
# Use with: tclsh drillBranchedPore.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

# Input:
set pdbIn block_hex.pdb
# Output:
set pdbOut branch.pdb

# Determine whether the position {x y z} is inside of pore and
# should be deleted.
proc insidePore {x y z} {
	set r0 8.
	set r1 5.
	
	if {$z < 0} {
		set isIn [expr $x*$x + $y*$y < $r0*$r0]
		return $isIn
	} else {
		# NOTE: Add your criteria to the expr commands.
		set isIn0 [expr ]
		set isIn1 [expr ]
		return [expr $isIn0 || $isIn1]
	}
}

proc drillPore {pdbIn pdbOut} {
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
		if {![insidePore $x $y $z]} {
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

drillPore $pdbIn $pdbOut



