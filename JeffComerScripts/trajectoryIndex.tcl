# List files associated with trajectories.
# Author: Jeff Comer <jcomer2@illinois.edu>

# Output:
set outFile pore2.0_index.txt

set count 0
set cd non_2ns_
set ad hairpin-E
set dd nonads_

set structDir(loop) /Projects/alek/jcomer/myhairpin/loop_first
set structName(loop) nw_pore+dna-all
set xscName(loop) run1_4V.restart
set dcdDir(loop) /Scr/nanopore/jcomer/myhairpin/loop_first_dcd

set structDir(coil) /Projects/alek/jcomer/myhairpin/nonstick_tail/1_run
set structName(coil) nw_pore+dna_E4V
set xscName(coil) steady_4V2.restart
set dcdDir(coil) /Scr/nanopore/jcomer/myhairpin/nonstick_tail_dcd

set structDir(coil_alek) /Projects/alek/jcomer/myhairpin/alek
set structName(coil_alek) nw_pore+dna-all
set xscName(coil_alek) hairpin-E4V-6.restart
set dcdDir(coil_alek) /Scr/nanopore/jcomer/myhairpin/alek_dcd

set structDir(double) /Projects/alek/jcomer/myhairpin/double
set structName(double) nw_pore+dna-all
set xscName(double) nonads_3V0.restart
set dcdDir(double) /Scr/nanopore/jcomer/myhairpin/double_dcd

set structDir(double_alek) /Projects/alek/jcomer/myhairpin/alek/dsDNA
set structName(double_alek) pore+dna-all
set xscName(double_alek) E5V-DS4-4.restart
set dcdDir(double_alek) /Scr/alek2/20_DS4/output/

proc expandRange {range} {
    set s [concat $range]
    set ret {}
    
    set last 0
    set inRange 0
    foreach n $s {
	if {$inRange && [string is integer $n]} {
	    for {set j [expr $last+1]} {$j <= $n} {incr j} {
		lappend ret $j
	    }
	    set inRange 0

	} elseif {[string is integer $n]} {
	    set last $n
	    lappend ret $n
	} elseif {[string equal $n to]} {
	    set inRange 1
	} else {
	    lappend ret $n
	    set inRange 0
	}
    }
    return $ret
}

proc extractVoltage {name} {
    set voltage [regexp -inline {\_[0-9.]*V} $name]
    return [string range $voltage 1 end]
}

proc add {class name dcdName dcdSet dcdFreq {offset 0.0}} {
    global out count structDir structName xscName dcdDir

    # Parse the dcdName, making substitutions.
    set n [string length $dcdName]
    set s ""
    for {set j 0} {$j < $n} {incr j} {
	# Get the current character.
	set c [string index $dcdName $j]
	
	if {[string equal $c "%"]} {
	    # Perform the substitutions.
	    set cNext [string index $dcdName [expr $j+1]]
	    switch $cNext {
		n {set s ${s}${name}}
		w {set s ${s}nw_}
		v {set s ${s}[extractVoltage ${name}]}
	    }
	    # Skip the next character.
	    incr j
	} else {
	    # Just add this character.
	    set s ${s}${c}
	}
    }
    
    set range [expandRange $dcdSet]
    
    puts "$count $name"
    puts $out "$count $class $name $structDir($class) $structName($class) $xscName($class) $dcdDir($class) $s \{$range\} $dcdFreq $offset"
    incr count
}

set out [open $outFile w]
puts $out "\# Trajectory analysis"
puts $out "\# index class name structDir structName xscName dcdDir dcdName dcdSet dcdFreq offset"

add loop loop_first_1V %w%n "0 to 2" 5000
add loop loop_first_1.5V %w%n "0 to 11" 5000
add loop loop_first_2V  %w%n "0 to 3" 5000 
add loop hairpin_first_1V  %w%n "0 to 10" 5000
add loop hairpin_first_1.5V  %w%n "0 to 3" 5000
add loop hairpin_first_2V  %w%n "0 to 6" 5000
add loop hairpin_first_2.5V  %w%n "0 to 6" 5000
add loop hairpin_first_3V %w%n "0 to 2" 5000
add loop hairpin_first_4V %w%n "0 to 2" 5000
add loop hairpin_first_6V %w%n "0 to 2" 5000
add coil coil_first_1V %w$cd%v 0 5000
add coil coil_cont_1V %w%n "0 to 8" 5000 9.0
add coil coil_first_1.5V %w$cd%v "1 to 4" 5000 3.36
add coil coil_cont_1.5V %w%n "0 to 6" 5000 10.3
add coil coil_first_2V %w$cd%v "0 to 7" 5000
add coil coil_first_3V %w$cd%v "0 to 2" 5000
add coil coil_first_4V %w$cd%va "{} 1 to 5" 5000 -2.7
add coil_alek coil_alek_4V %w$ad%v- "1 to 9" 5000
add coil_alek coil_alek_6V %w$ad%v- "1 to 5" 5000 -4.14
add double double_3V %w$dd%v "0 to 3" 5000
add double double_4V %w$dd%v "0 to 2" 5000
add double_alek double_alek_4V E4V-DS4-Runs1-5_50ps "" 50000
add double_alek double_alek_5V E5V-DS4-Runs1-4_50ps "" 50000
add double_alek double_alek_6V E6V-DS4-Runs1-4_50ps "" 50000
add double_alek double_alek_8V E8V-DS4-Runs1-3_50ps "" 50000
close $out



