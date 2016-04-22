# Make extrabond file for MGHH stability
# Usage: tclsh mghh_extrabonds.tcl pdbFile extrabondFile
# Author: Chen-Yu Li <cli56@illinois.edu> 
# 2015/5/19

set pattern {ATOM  (.{5}).(.{4}).(.{3}).(.)(.{4})....(.{8})(.{8})(.{8}).{6}.{6}.{6}(.{4})}

set count 0

set inFile [open [lindex $argv 0] r]
#set outFile [open [lindex $argv 1] w]
foreach line [split [read $inFile] \n] {

    if {[regexp $pattern $line matched index name resname chain resid x y z segname]} {

        incr count

        if {[regexp {MG} $name matched]} {

            puts [format "bond %8d %8d %8.3f %8.3f" [expr $count-1] $count 5000 1.94]
            puts [format "bond %8d %8d %8.3f %8.3f" [expr $count-1] [expr $count + 3]  5000 1.94]
            puts [format "bond %8d %8d %8.3f %8.3f" [expr $count-1] [expr $count + 6]  5000 1.94]
            puts [format "bond %8d %8d %8.3f %8.3f" [expr $count-1] [expr $count + 9]  5000 1.94]
            puts [format "bond %8d %8d %8.3f %8.3f" [expr $count-1] [expr $count + 12] 5000 1.94]
            puts [format "bond %8d %8d %8.3f %8.3f" [expr $count-1] [expr $count + 15] 5000 1.94]
        }


    }
}

close $inFile
#close $outFile


