# Compute the solvent-accessible surface area for a trajectory.
# to use: vmd -dispdev text -e sasaDcd.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

# requires {name dir struct dcdPrefix stride dcdSet}

proc compute {name dir struct dcdSet dcdFreq startFrame} {
    # Input:
    set psf $dir/${struct}.psf
    set pdb $dir/${struct}.pdb
    # Output:
    set outPrefix output/dmmp_orient_${name}
    # Parameters:
    set selText0 "resname DMMP and name P1"
    set selText1 "resname DMMP and name OD"
    set timestep 1.0
    set displayPeriod 200
    set stride 1

    set degrees [expr 180.0/(4.0*atan(1.0))]

    # Get the time change between frames in femtoseconds.
    set dt [expr $timestep*$dcdFreq*$stride]
    set dcdEnd [llength $dcdSet]

    # Open the output file.
    set out [open $outPrefix.dat w]

    # Load the system.
    mol load psf $psf pdb $pdb
    set sel0 [atomselect top $selText0]
    set sel1 [atomselect top $selText1]

    # Loop over the dcd files.
    set nFrames0 0
    foreach dcd $dcdSet {
	# Load the trajectory.
	animate delete all
	mol addfile $dcd first $startFrame step $stride waitfor all
	set nFrames [molinfo top get numframes]
	puts [format "Reading %i frames." $nFrames]

	for {set f 0} {$f < $nFrames} {incr f} {
	    molinfo top set frame $f

	    set n 0
	    set r0 [measure center $sel0 weight mass]
	    set r1 [measure center $sel1 weight mass]
	    set d [vecsub $r1 $r0]
	    set d [vecscale [expr 1.0/[veclength $d]] $d]

	    set z [lindex $r0 2]
	    #set a [expr $degrees*acos([lindex $d 2])]
	    set a [lindex $d 2]
	    
	    puts $out "$z $a"
	    incr n
	    
	    # Get the time in nanoseconds for this frame.
	    set t [expr ($nFrames0+$f)*$dt*1.e-6]
	    if {$f % $displayPeriod == 0} {puts "FRAME $f: $z $a"}
	}
	set nFrames0 [expr $nFrames+$nFrames0]
    }
    close $out

    $sel0 delete
    $sel1 delete
    mol delete top

    return
}
