# Make a new psf file with the atom records of one
# but the bonds, angles, and dihedrals of another.
# Useful for reconstituting psf files that have been
# edited with VMD.
# Use with: tclsh insertPsfAtomRecords.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

# Input:
set recordPsf cg_hairpin_built1_iterate.psf
set completePsf cg_hairpin_built1.psf
# Output
set outPsf cg_hairpin_iterate.psf

# Extract the atom records from the record psf.
proc extractAtomRecords {recordPsf} {
    set in [open $recordPsf r]
    set reading 0
    set recordList {}
    foreach line [split [read $in] \n] {
	if {$reading} {
	   # Quit if we have left the atoms section of the psf.
	   if {[string match "*!NBOND*" $line]} {break}
	   # Check that this is a valid line.
	   #if {[string length $line] < 10} {continue}

	   lappend recordList $line
       } elseif {[string match "*!NATOM*" $line]} {
	   # Don't store the data until we reach the atom records.
	   set reading 1
	   lappend recordList $line
       }

	# Skip to the atom records.
	if {[string match "*!NATOM*" $line]} {set reading 1}
    }
    close $in
    return $recordList
}

set recordList [extractAtomRecords $recordPsf]
puts "Read [llength $recordList] atom records from $recordPsf."

set out [open $outPsf w]
set in [open $completePsf r]
set insideRecords 0
foreach line [split [read $in] \n] {
    if {[string match "*!NATOM*" $line]} {
	# We are in the atom records.
	# Copy them in.
	set insideRecords 1
	puts "\nThese two lines usually should match: "
	puts "$recordPsf: [lindex $recordList 0]"
	puts "$completePsf: $line"

	foreach record $recordList {
	    puts $out $record
	}

    } elseif {[string match "*!NBOND*" $line]} {
	set insideRecords 0
    }

    # Write the line if we are not in the atom records.
    if {!$insideRecords} {
	puts $out $line
    }
}
close $in
close $out

puts "$outPsf was written successfully."
exit



