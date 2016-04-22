proc get_first_ts { xscfile } {
    set fd [open $xscfile r]
    gets $fd
    gets $fd
    gets $fd line
    set ts [lindex $line 0]
    close $fd
    return $ts
}
