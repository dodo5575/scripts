# Add harmonic constraints to silicon nitride.
# to use: vmd -dispdev text -e constrainSilicon.tcl
# jcomer2@uiuc.edu

# Parameters:
set sys at
set ion pot
set selText "serial 11980" 
set moveText "serial 11980"
set runCount 6
set occupancy 1.0
set stride 1
set coorFileList [glob /projects/jcomer/basepair/output/b*_at_*.coor]
# Input:
set indexFile diff_window.txt
set psf run_at.psf
set coor run_at.pdb
# Output:
set outPrefix ../init/init_${ion}_${sys}_pos

# Load the system.
mol load psf $psf
mol addfile $coor
set coorNum [llength $coorFileList]

# Just read a space delimited data file.
proc readData {fileName} {
    set in [open $fileName r]
    
    set r {}
    while {[gets $in line] >= 0} {
	if {[string match "#*" $line]} {continue}
	if {[string length $line] < 2} {continue}

	set tok [concat $line]
	set t [lindex $tok 0]
	lappend r $tok
    }

    close $in
    return $r
}


# Read the index file.
foreach zero {0} {set data [readData $indexFile]}
set nodePos {}
set nodeIndex {}
foreach d $data {
    lappend nodePos [lrange $d 1 3]
    lappend nodeIndex [lindex $d 0]
}

# Make the selections.
set selAll [atomselect top all]
set sel [atomselect top $selText]
set selMove [atomselect top $moveText]
puts "Performing WHAM on [$sel num] atoms."

#set closeFile ${dcdPrefix}[lindex $dcdSet 0]${dcdSuffix}
# Set the initial search conditions.
set n [llength $nodePos]

# Write the smd files.
foreach pos $nodePos i $nodeIndex {
    puts "System $i, position $pos"
    
    for {set run 0} {$run < $runCount} {incr run} {
	# Set the coordinates from a random coordinate file.
	animate delete all
	set coorInd [expr {int(floor(rand()*$coorNum))}]
	set coorFile [lindex $coorFileList $coorInd] 
	puts "USING coordinates from `$coorFile'."
	mol addfile $coorFile

	$selAll set occupancy 0.0
	$sel set occupancy $occupancy
	set currPos [measure center $sel weight mass]
	$selMove moveby [vecsub $pos $currPos]

	# Write the WHAM file.
	$selAll writepdb ${outPrefix}${i}-${run}.pdb
    }
}

$sel delete
$selAll delete
mol delete top
exit
