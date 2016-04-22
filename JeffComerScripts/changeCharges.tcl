# This script will change the types of atoms.
# Use with: vmd -dispdev text -e changeTypes.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set posCharge 0.9
set posType SI
set negCharge -0.45
set negType OSI
#Input:
set psf silica_trap.psf
set pdb silica_trap.pdb
#Output:
set finalPsf silica_trap_charges.psf
set tempPsf tmp.psf

mol load psf $psf pdb $pdb
set all [atomselect top all]

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

proc insertAtomRecords {recordList completePsf outPsf} {
    set out [open $outPsf w]
    set in [open $completePsf r]
    set insideRecords 0
    foreach line [split [read $in] \n] {
	if {[string match "*!NATOM*" $line]} {
	    # We are in the atom records.
	    # Copy them in.
	    set insideRecords 1
	    puts "\nThese two lines usually should match: "
	    puts "[lindex $recordList 0]"
	    puts "$line"

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
    return
}

# Set the charges.
set sel [atomselect top "type $posType"]
$sel set charge $posCharge
puts "Set the charges of [$sel num] atoms to $posCharge."
$sel delete
set sel [atomselect top "type $negType"]
$sel set charge $negCharge
puts "Set the charges [$sel num] atoms to $negCharge."
$sel delete

$all writepsf $tempPsf
$all delete

# Insert the atom records.
foreach zero {0} {set recordList [extractAtomRecords $tempPsf]}
insertAtomRecords $recordList $psf $finalPsf
puts "Wrote the complete psf $finalPsf."

exit



