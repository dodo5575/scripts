# Rename a set of segments.
# Use with: tclsh renameSegment.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set selText "all"
set subIn "U0 "
set subOut "U0_"
# Parameters:
set fileNamePrefix anneal_surface
# Input:
set psfIn ${fileNamePrefix}.psf
set pdbIn ${fileNamePrefix}_cen.pdb
# Output:
set psfOut ${fileNamePrefix}_better.psf
set pdbOut ${fileNamePrefix}_better.pdb

set letList "A B C D E F G H I J K L M N O P"
set segNameIn {}
set segNameOut {}
foreach let $letList {
    lappend segNameIn ${subIn}${let}
    lappend segNameOut ${subOut}${let}
}

proc renameSegmentPdb {segNameIn segNameOut pdbIn pdbOut} {
	# Open the pdb to extract the atom records.
	set count 0
	set out [open $pdbOut w]
	set in [open $pdbIn r]
		
	foreach line [split [read $in] \n] {
    		set string0 [string range $line 0 3]
		
		# Just write any line that isn't an atom record.
		if {![string match $string0 "ATOM"]} {
			puts $out $line
			continue
		}
		
		# Get the segment name.
		set segName [string trim [string range $line 72 75]]
		
		# Does this segment match?
		set segNameNew $segName
		foreach si $segNameIn so $segNameOut {
			if {[string equal $segName $si]} {
				set segNameNew [string range [format "%s    " $so] 0 3]
				break
			}
		}
		
		# Just write any atom record that doesn't
		# have the segment name we're changing.
		if {[string equal $segName $segNameNew]} {
			puts $out $line
			continue
		}
		
		# Generate the new pdb line.
		set temp0 [string range $line 0 71]
		set temp1 [string range $line 76 end]
				
		# Write the new pdb line.
		puts $out ${temp0}${segNameNew}${temp1}
		incr count
	}
	close $in
	close $out
	return $count
}

proc renameSegmentPsf {segNameIn segNameOut psfIn psfOut} {
	# Open the pdb to extract the atom records.
	set count 0
	set out [open $psfOut w]
	set in [open $psfIn r]
	
	set record 0
	set n 0
	set num 1
	
	foreach line [split [read $in] \n] {
    		# If we have finished with the atom records, just write the line.
		if {$n >= $num} {
			puts $out $line
			continue
		}
		
		if {!$record} {
			# Check if we have started the atom records.			
			if {[string match "*NATOM" $line]} {
				set record 1
				set numIndex [expr [string last "!" $line]-1]
				set num [string trim [string range $line 0 $numIndex]]
			}
			
			# Write the line.
			puts $out $line	
		} else {
			incr n			
			set segName [string trim [string range $line 9 12]]
			
			# Does this segment match?
			set segNameNew $segName
			foreach si $segNameIn so $segNameOut {
				if {[string equal $segName $si]} {
					set segNameNew [string range [format "%s    " $so] 0 3]
					break
				}
			}
		
			if {![string equal $segName $segNameNew]} {
				set temp0 [string range $line 0 8]
				set temp1 [string range $line 13 end]
				
				# Write the new line.
				puts $out ${temp0}${segNameNew}${temp1}
				incr count
			} else {
				# Just write the line.
				puts $out $line
			}
			
		}	
	}
	close $in
	close $out
	return $count
}

set nPsf [renameSegmentPsf $segNameIn $segNameOut $psfIn $psfOut]
puts "Changed $nPsf segment names in `$psfIn', producing `$psfOut'."
set nPdb [renameSegmentPdb $segNameIn $segNameOut $pdbIn $pdbOut]
puts "Changed $nPdb segment names in `$pdbIn', producing `$pdbOut'."
exit




