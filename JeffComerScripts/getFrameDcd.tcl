# Calculate the center of mass for a trajectory.
# Writes step time(ns) and position (nm) to a file.
# to use: vmd -dispdev text -e trackPositionDcd.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set pickFrame 500

# Input:
set psf ramp_ade16_f199.psf
set pdb ramp_ade16_f199.pdb
set dcdPrefix /Scr/nanopore/jcomer/single/single_dcd/ramp_ade16_
set dcdSuffix ".dcd"
set dcdSet {more_2ns0 more_2ns2 more_2ns3}
# Output:
set outPrefix ramp_ade_

# Load the system.
mol load psf $psf pdb $pdb
set sel [atomselect top all]

# Loop over the dcd files.
set nFrames0 0
set dcdEnd [llength $dcdSet]
for {set dcd 0} {$dcd < $dcdEnd} {incr dcd} {
    # Load the trajectory.
    animate delete all
    set dcdNum [lindex $dcdSet $dcd]
    mol addfile "${dcdPrefix}${dcdNum}${dcdSuffix}" type dcd waitfor all first $pickFrame last $pickFrame

    $sel writepdb ${outPrefix}${dcdNum}_f0.5ns.pdb
}
exit



