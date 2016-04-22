# Add harmonic constraints to silicon nitride.
# to use: vmd -dispdev text -e constrainSilicon.tcl
# jcomer2@uiuc.edu

# Parameters:
set selTextList {"index 20470" "index 20471" "index 20457" "index 20458"}
set geoList {{-1 0 0} {0 0 0} {0 1 0} {1 1 0}}
set occupancy 1.0
set stride 1
set maxShift 1000
# Input:
set indexFile win_dense.txt
set psf ../solution.psf
set coordPrefix scramble/scramble
# Output:
set outPrefix ../pdb/ions_pos

# Load the system.
mol load psf $psf
mol addfile ${coordPrefix}0.pdb

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
set selList {}
foreach st $selTextList {
    lappend selList [atomselect top $st]
}
puts "Performing WHAM on [llength $selList] selections."

#set closeFile ${dcdPrefix}[lindex $dcdSet 0]${dcdSuffix}
# Set the initial search conditions.
set n [llength $nodePos]

# Write the smd files.    
foreach pos $nodePos i $nodeIndex {
    puts "System $i, position $pos"

    animate delete all
    mol addfile ${coordPrefix}${i}.pdb
    
    $selAll set occupancy 0.0
    foreach s $selList geo $geoList {
	set r [vecscale [lindex $pos 2] $geo]

	set currPos [measure center $s weight mass]
	$s moveby [vecsub $r $currPos]
    }

    # Write the WHAM file.
    $selAll writepdb ${outPrefix}${i}.pdb
}

foreach s $selList {
    $s delete
}
$selAll delete
mol delete top
exit

