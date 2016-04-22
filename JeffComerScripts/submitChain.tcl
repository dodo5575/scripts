set inFile layer9.sh
set inp [open $inFile r]
set wallTime 12
set command "multinamdgo.sh"

set count 0
while {[gets $inp line] >= 0} {
    if {[string length $line] <= 1} {continue}
    if {[string match "\#*" $line]} {continue}
    
    if {$count == 0} {
	#puts "$cmd $line"
    } else {
	# Get the last line of the log file.
	set lastLog [concat [exec tail -n1 $env(HOME)/.multinamd.log]]
	set lastJob [lindex $lastLog 0]
	set cmd "exec $command -w 12 -d $lastJob $line"
	puts $cmd
	eval $cmd
    }

    #after 500
    incr count
}

close $inp
