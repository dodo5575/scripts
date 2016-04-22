# Add harmonic constraints to silicon nitride.
# to use: vmd -dispdev text -e constrainSilicon.tcl
# jcomer2@uiuc.edu

# Parameters:
set selText "name P1"
set moveText "resname DMMP"
set occupancy 1.0
set stride 1
set maxShift 1000
# Input:
set indexFile win_pos_inter.txt
set psf dmmp+mem.psf
set coor phos_pos0.pdb
# Output:
set outPrefix pdb1/phos_pos

# Load the system.
mol load psf $psf
mol addfile $coor

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
    
    $selAll set occupancy 0.0
    $sel set occupancy $occupancy
    set currPos [measure center $sel weight mass]
    $selMove moveby [vecsub $pos $currPos]

    # Write the WHAM file.
    $selAll writepdb ${outPrefix}${i}.pdb
}

$sel delete
$selAll delete
mol delete top
exit
