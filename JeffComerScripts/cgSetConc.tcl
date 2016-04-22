# Get the system size using measure minmax.
# to use: vmd -dispdev text -e getSize.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set conc 0.4
set margin 6.0
set posName SOD
set negName CL
set waterName H2O
#set qp 0.7; # Resetting the charges for now.
#set qn -0.7;
# Input:
set psf cg_hemo_ions1.0.psf
set pdb cg_hemo_ions1.0.pdb
# Output:
set finalPsf cg_hemo_${conc}M.psf
set finalPdb cg_hemo_${conc}M.pdb

set molar 55.523
set tempPsf tmp.psf
set otherText "not within $margin of (not name $waterName $posName $negName)"
set atomFormat [list name resname type mass charge element]

# Get the atom data given a name.
proc getAtomData {name} {
    global atomFormat
    set s [atomselect top "name $name"]
    set data [lindex [$s get $atomFormat] 0]
    $s delete 

    return $data
}

# Transform one kind of atom into another using
# atom data of the format from getAtomData.
proc replaceAtom {srcName destData otherText} {
    global atomFormat

    # Get the indices of atoms of $srcName.
    set sel [atomselect top "(name $srcName) and ($otherText)"]
    if {[$sel num] <= 0} {
	puts "Warning replaceAtom: No more valid $srcName atoms exist."
	return 0
    }
    set indexList [$sel get index]
    $sel delete
    set n [llength $indexList]

    # Choose an atom randomly.
    set j [expr int(floor(rand()*$n))]
    if {$j >= $n} {set j [expr $n-1]}

    # Copy the destination attributes into the source atom.
    set index [lindex $indexList $j]
    set s [atomselect top "index $index"]
    $s set $atomFormat [list $destData]
    $s delete

    return 1
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
    return
}

# Get the water data.
mol load psf $psf pdb $pdb
set waterSel [atomselect top "name $waterName"]
set waterCount [$waterSel num]
set waterN [expr 4*$waterCount]
set waterData [getAtomData $waterName]
puts "Number of water beads: $waterCount"
puts "Number of waters (hypothetically): [expr $waterN]"
$waterSel delete

# Select the ions.
set posSel [atomselect top "name $posName"]
#$posSel set charge $qp
set posN [$posSel num]
set posQ [lindex [$posSel get charge] 0]
set posData [getAtomData $posName]
$posSel delete
set negSel [atomselect top "name $negName"]
#$negSel set charge $qn
set negN [$negSel num]
set negQ [lindex [$negSel get charge] 0]
set negData [getAtomData $negName]
$negSel delete

# Get the charge that must be accounted for.
set otherSel [atomselect top "not name $posName $negName"]
set otherCharge [measure sumweights $otherSel weight charge]
$otherSel delete

# Get the number of counterions.
if {$otherCharge < 0} {
    set posCtr [expr int(floor(-$otherCharge/$posQ))]
    set negCtr 0
} else {
    set negCtr [expr int(floor(- $otherCharge/$negQ))]
    set posCtr 0
}

# Get the number necessary for the desired concentration.
set c [expr $conc/$molar]
set waterNum [expr ($waterN + 4*($posN + $negN - $posCtr - $negCtr))/(1+8*$c)]
set posNum [expr int(floor($waterNum*$c + $posCtr))]
set negNum [expr int(floor($waterNum*$c + $negCtr))]
puts "\nInitial $waterName number: $waterN" 

puts "Initial $posName number: $posN" 
puts "Initial $negName number: $negN"
puts "\nFinal $waterName number: $waterNum" 
puts "Final $posName number: $posNum" 
puts "Final $negName number: $negNum" 

puts "Replacing atoms..."
# Replace ions or water to get desired concentration.
if {$posNum > $posN} {
    # There are too few ions.
    set n [expr $posNum - $posN]
    puts "Replacing $n $waterName atoms with $posName atoms..."
    for {set i 0} {$i < $n} {incr i} {
	replaceAtom $waterName $posData $otherText
    }
} else {
    # There are too many ions.
    set n [expr $posN - $posNum]
    puts "Replacing $n $posName atoms with $waterName atoms..."
    for {set i 0} {$i < $n} {incr i} {
	replaceAtom $posName $waterData $otherText
    }
}

# Replace ions or water to get desired concentration.
if {$negNum > $negN} {
    # There are too few ions.
    set n [expr $negNum - $negN]
    puts "Replacing $n $waterName atoms with $negName atoms..."
    for {set i 0} {$i < $n} {incr i} {
	replaceAtom $waterName $negData $otherText
    }
} else {
    # There are too many ions.
    set n [expr $negN - $negNum]
    puts "Replacing $n $negName atoms with $waterName atoms..."
    for {set i 0} {$i < $n} {incr i} {
	replaceAtom $negName $waterData $otherText
    }
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

# Check the charge.
set q [measure sumweights $all weight charge]
puts "\nTotal charge: $q"
$all delete

# Check the concentration.
set selTextList [list "name $posName" "name $negName"]
set waterSel [atomselect top "name $waterName"]
set waterNumber [expr 4*[$waterSel num]]
$waterSel delete
foreach s $selTextList {
    set sel [atomselect top $s]
    set n [$sel num]
    puts "Number of $s: $n"
    puts "Concentration of $s: [expr $molar*$n/$waterNumber] mol/kg"
    $sel delete
}

mol delete top
exit



