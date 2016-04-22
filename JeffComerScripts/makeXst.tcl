# Make a NAMD xst file by replicating a NAMD xsc file.
# Use with: tclsh makeXst.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set stride 25000
set n 1000
# Input:
set xsc hairpin-E4V-6.restart.xsc
# Output:
set xst replicate.xst

# Get the configuration.
set in [open $xsc r]
foreach line [split [read $in] \n] {
	if {[string match "#*" $line]} continue;
	if {[string length $line] < 2} continue;
	set config [split $line]
}
close $in

# Write the trajectory.
set out [open $xst w]
puts $out "# NAMD extended system configuration restart file"
puts $out "#\$LABELS step a_x a_y a_z b_x b_y b_z c_x c_y c_z o_x o_y o_z"
for {set i 1} {$i <= $n} {incr i} {
	puts -nonewline $out "[expr $i*$stride]"
	for {set j 1} {$j < 11} {incr j} {
		puts -nonewline $out " [lindex $config $j]"
	}
	puts $out ""
}
close $out





