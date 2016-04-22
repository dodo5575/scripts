# This script will remove all residues in the selection
# from psf and pdf files.
# Use with: vmd -dispdev text -e removeResidues.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set conc 1.0
set waterText "name H2O"
set posData {"QD" "SOD" "SOD" "ION" 0.7}
set negData {"QA" "CL" "CL" "ION" -0.7}
#Input:
set psf water_box.psf
set pdb water_box_cen.pdb
#Output:
set finalPsf cg_ions_${conc}M.psf
set finalPdb cg_ions_${conc}M.pdb

set molar 55.523
set margin 2.0
set tempPsf tmp.psf

proc addIon {ionData waterText margin} {
    set sel [atomselect top "($waterText) and not within $margin of not ($waterText)"]
    set indexList [$sel get index]
    $sel delete
    set n [llength $indexList]

    set j [expr int(floor(rand()*$n))]
    if {$j >= $n} {set j [expr $n-1]}

    set index [lindex $indexList $j]
    set s [atomselect top "index $index"]
    $s set {type name resname segname charge} [list $ionData]
    $s delete
}

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
}

mol load psf $psf pdb $pdb
set water [atomselect top $waterText]
set waterN [$water num]
set waterCount [expr 4*$waterN]
set c [expr $conc/$molar]

set n [expr int(floor($c/(1.0+8.0*$c)*$waterCount))]
puts "mole ratio: $c"
puts "number of water beads: $waterN"
puts "number of waters: $waterCount"
puts "final number of waters: [expr $waterCount - 2*$n]"
puts "number of ions: $n"

# Add the ions.
for {set i 0} {$i < $n} {incr i} {
    addIon $posData $waterText $margin
    addIon $negData $waterText $margin
}

# Write the result.
set all [atomselect top all]
$all writepsf $tempPsf
$all writepdb $finalPdb
puts "Wrote $finalPdb."

# Insert the atom records.
foreach zero {0} {set recordList [extractAtomRecords $tempPsf]}
insertAtomRecords $recordList $psf $finalPsf
puts "Wrote $finalPsf."

$all delete
mol delete top
exit



