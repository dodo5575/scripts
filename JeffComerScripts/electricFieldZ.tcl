# Compute the electric field along the z-axis from a
# ".dx" file containing the potential.
# The z-basis vector must be orthogonal to the other two.
# Use with: tclsh electricFieldZ.tcl input_dx
# Author: Jeff Comer <jcomer2@illinois.edu>

if {$argc != 1} {
	puts "This script requires one argument: "
	puts "  input .dx file specifying the potential"
	puts ""
	exit
}

# Input:
set inData [lindex $argv 0]
# Output:
set outTmp ${inData}.tmp
set outData ${inData}.Ez

# Compute the negative derivative along z for
# a list a values with z fast, y medium, and x slow.
proc diff {potentialVar dz lx ly lz} {
	upvar $potentialVar potential
	
	set eFieldZ {}
	set lz1 [expr $lz - 1]
	set n 0
	set k 0
		
	# Get the neighbor for the first node.
	set pot0 [lindex $potential $lz1]
	
	foreach pot $potential {
		lappend eFieldZ [expr ($pot0 - $pot)/$dz]
		set pot0 $pot
		
		# Compensate for periodic boundaries.
		set n [expr $n + 1]
		set k [expr $k + 1]
		if {$k >= $lz} {
			set k 0
			set pot0 [lindex $potential [expr $n + $lz1]]
		}
	}
	
	return $eFieldZ
}

# Write a list of variables in .dx array format.
proc writeArray {data outStream} {
	foreach {a b c} $data {
		puts $outStream [string trimright "$a $b $c"]
	}
}

# Shift the origin of the output file to obtain O(dz^2) accuracy.
proc shiftOriginZ {shift inName outName} {
	set in [open $inName r]
	set out [open $outName w]
	
	foreach line [split [read $in] \n] {
		if {[string match "origin *" $line]} {
			set tokens [split $line " "]
			foreach {word ox oy oz} $tokens {break}
			set oz [expr $oz + $shift]
			puts $out "origin $ox $oy $oz"
		} else {
			# Just write all lines not containing the origin.
			puts $out $line
		}
	}
	close $out
	close $in
}

# Write the header.
set in [open $inData r]
set out [open $outTmp w]
puts "Extracting potential data..."

# Get the data.
set data {}
set lx 0
set ly 0
set lz 0
set readingArray 0
set potential {}
set nItems 0
set n 0
foreach line [split [read $in] \n] {
	if {[string match "#*" $line]} {
		puts $out $line
	} elseif {$readingArray} {
		set tokens [split $line " "]
		# Add the potential values to the list.
		foreach t $tokens {
			lappend potential $t
		}
		
		set n [expr $n + [llength $tokens]]
		if {$n >= $nItems} {
			puts "Finished reading potential values.\n"
			set readingArray 0
			
			puts -nonewline "Calculating the electric field..."
			set eField [diff potential $dz $lx $ly $lz]
			puts "Done."
			puts -nonewline "Writing the electric field..."
			writeArray $eField $out
			puts "Done."
		}
	} elseif {[string match "delta *" $line]} {
		puts $out $line
		set tokens [split $line " "]
		set deltaZ [lindex $tokens 3]
		if {$deltaZ != 0} {
			set dz $deltaZ
			puts "Grid spacing along z is ${dz}."
			puts "For optimal accuracy we'll shift origin along z by [expr -0.5*$dz]."
		}
	} elseif {[string match "*gridpositions*" $line]} {
		puts $out $line
		set tokens [split $line " "]
		set lx [lindex $tokens 5]
		set ly [lindex $tokens 6]
		set lz [lindex $tokens 7]
		puts "The grid size is ($lx $ly $lz)."
	} elseif {[string match "*class array type double*" $line]} {
		puts $out $line
		set tokens [split $line " "]
		set nItems [lindex $tokens 9]
		if {$n == [expr $lx*$ly*$lz]} {
			puts "Warning: Grid dimensions do not agree with data size."
		}
		puts "Reading $nItems potential values."
		
		set readingArray 1
		set n 0
	} else {
		puts $out $line
	}
}
close $out
close $in
puts "The file ${outTmp} was written successfully.\n"

puts -nonewline "Shifting the origin..."
shiftOriginZ [expr -0.5*$dz] $outTmp $outData
puts "Done."
puts "The file ${outData} was written successfully."

exit



