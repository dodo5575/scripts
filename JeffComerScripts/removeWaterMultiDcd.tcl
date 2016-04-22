# Remove the water from dcd files.
# Author: Jeff Comer <jcomer2@illinois.edu>

set stride 1
set selText "(not water) and (not resname SIN)"
# Input:
set psf0 pore1.2_all.psf
set pdb0 pore1.2_all.pdb
set dcdDir /Scr/nanopore/jcomer/myhairpin/pore1.2_coil_first_dcd
set dcdPrefix coil_p1.2_4Va11
set dcdSet {.1 .3 .5}
set dcdSuffix ".dcd"
# Output:
set outPrefix nw_
set psf ${outPrefix}${psf0}
set pdb ${outPrefix}${pdb0}

# Load the structure and make the selection.
mol load psf $psf0 pdb $pdb0
set sel [atomselect top $selText]
set nAtoms [$sel num]
puts "\nNOTE: Retaining $nAtoms atoms defined by $selText."
animate delete all

set dcdEnd [llength $dcdSet]
for {set dcd 0} {$dcd < $dcdEnd} {incr dcd} {
    # Load the trajectory.
    set dcdNum [lindex $dcdSet $dcd]
    mol addfile "$dcdDir/${dcdPrefix}${dcdNum}${dcdSuffix}" type dcd step $stride waitfor all
    set nFrames [molinfo top get numframes]
    puts [format "Reading %i frames." $nFrames]
    set last [expr $nFrames - 1]
}


# Write the psf once.
$sel writepsf $psf
$sel writepdb $pdb

animate write dcd "$dcdDir/${outPrefix}${dcdPrefix}${dcdSuffix}" beg 0 end $last waitfor all sel $sel top

$sel delete
mol delete top
exit



