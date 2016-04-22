set stepFile pull_4V3.restart.xsc

# Restart using the timestep given in the xsc file.
proc getFirstStep { xscFile } {
	set in [open $xscFile r]
	foreach line [split [read $in] "\n"] {
		if {![string match "#*" $line]} {
	  		set param [split $line]
			if {[llength $param] > 0} {
				set ts [lindex $param 0]
				break
			}
		}
	}
	close $in
	return $ts
}

if [file exists $stepFile] {
	puts [getFirstStep $stepFile]
}



