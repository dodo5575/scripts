# Remove the water from dcd files.
# Author: Jeff Comer <jcomer2@illinois.edu>

set stride 1
set selText "(not water) and (not resname SIN)"
# Input:
set psf0 ecoRI_2V.psf
set pdb0 ecoRI_2V.pdb
#set dcdDir /Scr/nanopore/jcomer/ecoRI
set dcdDir /Scr/7day-2/jcomer
set dcdPrefix eco_grisha_2Va
set dcdSet {0 1 2 3}
set dcdSuffix ".dcd"
# Output:
set outDir /Scr/nanopore/jcomer/ecoRI
set outPrefix nw_
set psf ${outPrefix}${psf0}
set pdb ${outPrefix}${pdb0}

set window 400
# Load the structure and make the selection.
mol load psf $psf0 pdb $pdb0
set sel [atomselect top $selText]
set nAtoms [$sel num]
puts "\nNOTE: Retaining $nAtoms atoms defined by $selText."

# Write the psf.
animate delete all
set dcdNum [lindex $dcdSet 0]
mol addfile "$dcdDir/${dcdPrefix}${dcdNum}${dcdSuffix}" type dcd first 0 last 0 waitfor all
$sel writepsf $psf
$sel writepdb $pdb

set dcdEnd [llength $dcdSet]
for {set dcd 0} {$dcd < $dcdEnd} {incr dcd} {
    # Load the trajectory in chunks of length $window.
    set nFrames 1
    set dcdNum [lindex $dcdSet $dcd]
    for {set j 0} {$nFrames > 0} {incr j} {
	animate delete all
	# Load this part of the trajectory.
	set first [expr $j*$window]
	set last [expr $first+$window-1]
        mol addfile "$dcdDir/${dcdPrefix}${dcdNum}${dcdSuffix}" type dcd step $stride first $first last $last waitfor all

	# Write it to a file.
	set nFrames [molinfo top get numframes]
	if {$nFrames <= 0} {continue}
	puts [format "Reading %i frames." $nFrames]
	animate write dcd "$outDir/${outPrefix}${dcdPrefix}${dcdNum}.${j}${dcdSuffix}" beg 0 end [expr $nFrames-1] waitfor all sel $sel top
    }
}

$sel delete
mol delete top
exit



