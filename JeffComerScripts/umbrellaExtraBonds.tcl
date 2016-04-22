# Add harmonic constraints to silicon nitride.
# to use: vmd -dispdev text -e constrainSilicon.tcl
# jcomer2@uiuc.edu

# Parameters:
set nameList {pot-pot pot-chl chl-chl}
set pairList {{20470 20471} {20471 20457} {20457 20458}}
set occupancy 1.0
set stride 1
set maxShift 1000
# Input:
set indexFile win_dense.txt
# Output:


# Just read a space delimited data file.
proc readData {fileName} {
    set in [open $fileName r]
    
    set r {}
    while {[gets $in line] >= 0} {
	if {[string match "#*" $line]} {continue}
	if {[string length $line] < 2} {continue}

	set tok [concat $line]
	lappend r [lrange $tok 0 1]
    }

    close $in
    return $r
}

# Read the index file.
set data [readData $indexFile]
set nodePos {}
set nodeSpring {}
set nodeIndex {}
foreach d $data {
    lappend nodePos [lrange $d 1 3]
    lappend nodeSpring [lrange $d 4 6]
    lappend nodeIndex [lindex $d 0]
}

foreach pos $nodePos spring $nodeSpring index $nodeIndex {
    foreach name $nameList pair $pairList {
	set outFile ${name}_pos${index}.txt
	set out [open $outFile w]

	set ind0 [lindex $pair 0]
	set ind1 [lindex $pair 1]
	set k [expr 0.5*[lindex $spring 2]]
	set d [lindex $pos 2]

	puts $out "bond $ind0 $ind1 $k $d"
	close $out
    }
}