# Make NAMD configuration files.
# Use with: tclsch removeWater.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set prefix "uz"
set name eq0
set n 3

for {set i 0} {$i < $n} {incr i}  {
	set out [open "${prefix}${i}_${name}.csh" w]

	puts $out "#!/bin/csh"

	puts $out "#$ -e log-local-err.log"
	puts $out "#$ -o log-local-out.log"
	puts $out "#$ -cwd"
	
	puts $out "namd2 ${prefix}${i}_${name}.namd >! ${prefix}${i}_${name}.log"
	puts $out "exit 0"

	close $out
}




