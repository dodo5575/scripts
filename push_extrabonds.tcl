# Make extrabond file for MGHH stability
# Usage: vmd -dispdev text -e push_extrabonds.tcl -args ndxFile pdbFile extrabondFile
# Author: Chen-Yu Li <cli56@illinois.edu> 
# 2015/6/4

set pattern {ri ([0-9]*) ([0-9]*) & a P}

set pairs ""

set inFile [open [lindex $argv 0] r]
foreach line [split [read $inFile] \n] {

    if {[regexp $pattern $line matched r1 r2]} {

        set pair [list $r1 $r2]

        lappend pairs $pair

    }
}

close $inFile

mol new [lindex $argv 1]

set outFile [open [lindex $argv 2] w]
foreach pair $pairs {

    # VMD and NAMD are 0-based
    set r1 [expr [lindex $pair 0] -1]
    set r2 [expr [lindex $pair 1] -1]
    set index1 [[atomselect top "residue $r1 and name P"] get index]
    set index2 [[atomselect top "residue $r2 and name P"] get index]
    if {$index1 != "" && $index2 != ""} {

        puts $outFile [format "bond %8d %8d %8d %8d" $index1 $index2 1 30]        
    }
    
}
close $outFile

exit

