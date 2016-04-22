# Change a type
# Author: Jeff Comer <jcomer2@illinois.edu>

if {$argc != 3} {
	puts "This script requires three arguments: "
	puts "  input type"
	puts "  output type"
	puts "  input file"
	puts ""
	exit
}

set inType [lindex $argv 0]
set outType [lindex $argv 1]
# Input:
set inData [lindex $argv 2]
# Output:
set outData ${inData}.type

if {[file exists $outData]} {
    puts "File $outData exists. Overwrite?"
    gets stdin q
    if {[string match -nocase "y*" $q]} {
	puts "Overwriting $outData."
    } else {
	exit
    }
}

set headerList {BONDS THETAS PHI IMPHI NONBONDED END}
puts "\nScanning lines..."
set in [open $inData r]
set out [open $outData w]
set count 0
set cont 0
foreach line [split [read $in] \n] {
    # Ignore comments.
    if {[string match "!*" $line]} {continue}
    if {[string equal [string index $line 0] "*"]} {
	puts $line
	continue
    }

    # Write continuations.
    if {$cont} {
	puts $out $line
	set cont [string equal [string index $line end] "-"]
	continue
    }

    # Maintain breaks.
    if {[string length $line] < 1} {
	puts $out ""
	continue
    }
    
    # Just write the section headers.
    foreach header $headerList {
	if {[string match "${header}*" $line]} {
	    puts $out $line
	    set cont [string equal [string index $line end] "-"]
	    continue
	}
    }

    # Write the replacement lines.
    set n [regexp -all $inType $line]
    if {$n >= 0} {
	# Replace each instance of inType with outType one at a time.
	# Hence, we assume not interaction between outType atoms.
	set start 0
	for {set i 0} {$i < $n} {incr i} {
	    regexp -start $start -indices $inType $line ind
	    puts $out [regsub -start $start $inType $line $outType]
	    set start [expr [lindex $ind 1]+1]
	}

	incr count
    }
}

close $out
close $in

puts "Extracted $count lines."
puts "Wrote `$outData'."



