# Match one selection with another.
# Use with: vmd -dispdev text -e fit.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set azimuth -15.0
set latitude 20.0
set displace 1.5
set extraAzimuth 90
# Input:
set inPrefix ../1_build_all/struct/helix3_
# Output:
set outPrefix struct/inter3_

set pi [expr {4.0*atan(1.0)}]
#set theta [expr {$azimuth*$pi/180.0}]
#set phi [expr {$latitude*$pi/180.0}]

set seqList {}
foreach a {A T G C} {
    foreach b {A T G C} {
	foreach c {A T G C} {
	    lappend seqList ${a}${b}${c}
	}
    }
}

foreach seq $seqList {
    set inFile ${inPrefix}${seq}.pdb
    set outFile ${outPrefix}${seq}.pdb
    
    mol load pdb $inFile
   
    set all [atomselect top all]
    set resList [lsort -unique -integer [$all get resid]]

    set resA [lrange $resList 0 2]
    set resB [lrange $resList 3 5]

    for {set i 0} {$i < 3} {incr i} {
	set sign [expr {($i-1)}]
	
	set a [lindex $resA $i]
	set b [lindex $resB [expr {2-$i}]]
	set sel [atomselect top "resid $a $b"]
	$sel move [transaxis z [expr {$sign*$azimuth+$extraAzimuth}]]

	set cen [measure center $sel weight mass]
	$sel moveby [vecinvert $cen]
	$sel move [transaxis x $latitude]
	$sel moveby $cen
	$sel moveby [list 0 0 [expr {$displace*$sign}]]
	$sel delete
    }
    
    $all writepdb $outFile
    $all delete
}

exit
