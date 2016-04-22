# Calculate the number of charge carriers and current (nA) with multiple
# selections for a trajectory.
# Output format is t(ns) nSel0 ISel0(nA) nSel1 ISel1(nA) ... nTotal ITotal(nA).
# to use: vmd -dispdev text -e multicurrentDcd.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>
set cutoff 10.0
set dcdFreq 5000
set wallText "resname SIN"
set backboneHydrogen "H1' H2' H2'' H3' H4' H5' H5''"
set selName [list "backbone" "base"]
set selText [list "segname HAIR and name $backboneHydrogen" "segname HAIR and name \"H.*\" and not name $backboneHydrogen"] 
set timestep 1.0
# Input:
set pdb pore+dna-all.pdb
set psf pore+dna-all.psf
set dcdPrefix /Projects/alek/HAIRPIN/20_HAIR-PORE/5_RUNS/hairpin-
set dcdSuffix ".dcd"
set dcdSet {E6V-4 E6V-5}
set coords pore2.0_coords.txt
# Output:
set outFilePrefix stick_tail

# Do the cell decomposition.
set once {0}
foreach i $once {
source cellDecomposition.tcl
set pos [loadPoints $coords]
set n [llength $pos]
puts "\nCell decomposition"
puts "Performing a cell decomposition on $n points with a cutoff of $cutoff"
set gridDef [cellGrid $pos $cutoff]
set origin [lindex $gridDef 0]
set grid [lindex $gridDef 1]
puts "Origin: $origin"
puts "Grid: $grid"
set nCells [expr [lindex $grid 0]*[lindex $grid 1]*[lindex $grid 2]]
puts "Number of cells: $nCells"
set neigh [cellNeighbors $grid]
puts "Decomposing..."
cellDecompose $pos $cutoff $origin $grid cell
puts "Cell decomposition complete"
}

# Open the output files.
set out {}
foreach s $selName {
    lappend out [open ${outFilePrefix}_${s}.dat w]
}

set dcdEnd [llength $dcdSet]
# Get the time change between frames in femtoseconds.
set dt [expr $timestep*$dcdFreq]

# Load the system.
mol load psf $psf pdb $pdb
set sel {}
foreach s $selName st $selText {
    set ss [atomselect top $st]
    lappend sel $ss
    puts "$s: [$ss num] atoms"
}


# Loop over the dcd files.
set nFrames0 1
for {set dcd 0} {$dcd < $dcdEnd} {incr dcd} {

# Load the trajectory.
animate delete all
set dcdNum [lindex $dcdSet $dcd]
mol addfile "${dcdPrefix}${dcdNum}${dcdSuffix}" waitfor all
set nFrames [molinfo top get numframes]
puts [format "Reading %i frames." $nFrames]

# Move forward, computing
# current at each step.
for {set f 1} {$f < $nFrames} {incr f} {
    molinfo top set frame $f

    # Get the time in nanoseconds for this frame.
    set t [expr ($nFrames0+$f)*$dt*1.e-6]

    # Write the time.
    puts "FRAME $f: $t"

    foreach s $sel o $out {
	set selPos [$s get {x y z}]
	foreach r $selPos {
	    set index [cellLookup $r $cutoff $origin $grid]
	    set distance $cutoff
	    
	    # Find the nearest distance.
	    foreach n [lindex $neigh $index] {
		set points $cell($n)
		
		foreach p $points {
		    set d [veclength [vecsub $r $p]]
		    if {$d < $distance} {
			set distance $d
		    }
		}
	    }
	    if {$distance < $cutoff} {
		puts $o $distance
	    }
	}
    }
}
set nFrames0 [expr $nFrames+$nFrames0]
}

foreach o $out {
    close $o
}
exit




