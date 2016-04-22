# Based on continuity of resid, separate pdb into a list of segments
# Usage: tclsh pdb2chainsep.tcl pdbFile
# Author: Chen-Yu Li <cli56@illinois.edu> 
# 2015/5/6


set prefix "dna_chain"

set file_index 0

set resid_current 1

set chain "A"

set inStream [open [lindex $argv 0] r]
set out_current [open ${prefix}_${file_index}.pdb w]
foreach line [split [read $inStream] \n] {

    set string1 [string range $line 0 3]
    #ATOM

    if {[string equal {ATOM} $string1]} {
        set string4 [string trim [string range $line 17 19]]
        #resname

        set string5 [string trim [string range $line 22 25]]
        #resid

        set line [string replace $line 21 21 $chain]

        if {$string5 >= $resid_current} {

            set resid_current $string5
            puts $out_current $line


        } else {

            close $out_current
            puts [format "chain=%s" $chain]

            set resid_current 1

            incr file_index

            set chain [format %c [expr [scan $chain %c] + 1]]

            set out_current [open ${prefix}_${file_index}.pdb w]

        }
    }

}
puts [format "chain=%s" $chain]

close $inStream
close $out_current


