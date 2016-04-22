# Remove the water from dcd files.
# Author: Jeff Comer <jcomer2@illinois.edu>

set stride 1
set selText "(not water) and (not resname SIN)"
# Input:
set psf0 wad_all.psf
set pdb0 wad_all.pdb
set dcdDir /Scr/7day-2/jcomer
set dcdPrefix ""
set dcdSet {helix_0V0_scaling parmbsc0_0V0 parmbsc0_0V0_scaling}
set dcdSuffix ".dcd"
# Output:
set outDir /Scr/nanopore/jcomer/myhairpin/race1_dcd
set outPrefix nw_
set psf ${outPrefix}${psf0}
set pdb ${outPrefix}${pdb0}

# Load the structure and make the selection.
mol load psf $psf0 pdb $pdb0
set sel [atomselect top $selText]
set nAtoms [$sel num]
puts "\nNOTE: Retaining $nAtoms atoms defined by $selText."

set dcdEnd [llength $dcdSet]
for {set dcd 0} {$dcd < $dcdEnd} {incr dcd} {
    # Load the trajectory.
    animate delete all
    set dcdNum [lindex $dcdSet $dcd]
    mol addfile "$dcdDir/${dcdPrefix}${dcdNum}${dcdSuffix}" type dcd step $stride waitfor all
    set nFrames [molinfo top get numframes]
    puts [format "Reading %i frames." $nFrames]
    set last [expr $nFrames - 1]

    # Write the psf once.
    if {$dcd == 0} {
	$sel writepsf $psf
	$sel writepdb $pdb
    }

    animate write dcd "$outDir/${outPrefix}${dcdPrefix}${dcdNum}${dcdSuffix}" beg 0 end $last waitfor all sel $sel top
}

$sel delete
mol delete top
exit



