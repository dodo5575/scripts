# Author: Jeff Comer <jcomer2@illinois.edu>
set stride 200
set prefix traj_phant
set coorList [glob $prefix*]

# Extract the numbers so that we can sort by them.
set sortList {}
foreach coor $coorList {
    set num [string range [regexp -inline -all "\[.\]\[0-9\]*\[.\]" $coor] 1 end-1]
    lappend sortList [list $coor $num]
}

# Sort.
set sortList [lsort -integer -index 1 $sortList]

# Add the files.
foreach pair $sortList {
    set fileName [lindex $pair 0]
    mol addfile $fileName
}
