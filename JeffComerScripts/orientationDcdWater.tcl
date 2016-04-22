# Compute the solvent-accessible surface area for a trajectory.
# to use: vmd -dispdev text -e sasaDcd.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

# requires {name dir struct dcdPrefix stride dcdSet}

# Input:
set psf $dir/${struct}.psf
set pdb $dir/${struct}.pdb
# Output:
set outPrefix output/water_orient_${name}
# Parameters:
set selText "water"
set dcdFreq 5000
set timestep 1.0

set degrees [expr 180.0/(4.0*atan(1.0))]

# Get the time change between frames in femtoseconds.
set dt [expr $timestep*$dcdFreq]
set dcdEnd [llength $dcdSet]

# Open the output file.
set out [open $outPrefix.dat w]

# Load the system.
mol load psf $psf pdb $pdb
set sel [atomselect top $selText]
foreach zero {0} {set resList [lsort -unique [$sel get {segname resid}]]}

set selList0 {}
set selList1 {}
foreach res $resList {
    set s [lindex $res 0]
    set r [lindex $res 1]
    lappend selList0 [atomselect top "segname $s and resid $r and name OH2"]
    lappend selList1 [atomselect top "segname $s and resid $r and name H1 H2"]
}

# Loop over the dcd files.
set nFrames0 0
foreach dcd $dcdSet {
    # Load the trajectory.
    animate delete all
    mol addfile $dcd first $startFrame waitfor all
    set nFrames [molinfo top get numframes]
    puts [format "Reading %i frames." $nFrames]

    for {set f 0} {$f < $nFrames} {incr f $stride} {
	molinfo top set frame $f

	set n 0
	foreach s0 $selList0 s1 $selList1 {
	    set r0 [measure center $s0 weight mass]
	    set r1 [measure center $s1 weight mass]
	    set d [vecsub $r1 $r0]
	    set d [vecscale [expr 1.0/[veclength $d]] $d]

	    set z [lindex $r0 2]
	    #set a [expr $degrees*acos([lindex $d 2])]
	    set a [lindex $d 2]
	    
	    puts $out "$z $a"
	    incr n
	}
	
	# Get the time in nanoseconds for this frame.
	set t [expr ($nFrames0+$f)*$dt*1.e-6]
	if {$f % $displayPeriod == 0} {puts "FRAME $f: $z $a"}
    }
    set nFrames0 [expr $nFrames+$nFrames0]
}
close $out
$sel delete

foreach s0 $selList0 s1 $selList1 {
    $s0 delete
    $s1 delete
}
mol delete top
