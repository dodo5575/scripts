# Calculate the electric force from a grid for a trajectory.
# Writes step time(ns) and force (pN) to a file.
# to use: vmd -dispdev text -e electrostaticForce.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set samplesZ 400
set zStart -200
set zEnd 200
# Input:
set forceFile ez_open_p1.6_4V_center.dat
switch 1 {
    0 {
	set selText "segname HAIR"
	set psf hairpin_E4V.psf
	set coor coil_first_E4V.pdb
	set outputFile fc_coil_p1.6.dat
    }
    1 {
	set selText "segname HAIR"
	set psf hairpin_E4V.psf
	set coor loop_first_E4V.pdb
	set outputFile fc_loop_p1.6.dat
    }
    2 {
	set selText "segname ADNA BDNA"
	set psf double.psf
	set coor double_E4V.pdb
	set outputFile fc_double_p1.6.dat
    }
}

proc readForce {fileName} {
    set in [open $fileName r]
    set z {}
    set f {}

    foreach line [split [read $in] "\n"] {
	set tok [concat $line]
	
	# Ignore blank lines.
	if {[llength $tok] < 1} {continue}
	
	lappend z [lindex $tok 0]
	lappend f [lindex $tok 1]
    }
    return [list $z $f]
}

foreach null {0} {set data [readForce $forceFile]}
set dz [expr [lindex $data 0 1]-[lindex $data 0 0]]
set z0 [lindex $data 0 0]
set z1 [lindex $data 0 end]
foreach null {0} {set eFieldZ [lindex $data 1]}
set nz [llength $eFieldZ]

puts "Read $forceFile"
puts "Found $nz elements from $z0 to $z1 by $dz"

# Load the system.
mol load psf $psf
animate delete all
mol addfile $coor type pdb waitfor all
set sel [atomselect top $selText]
foreach null {0} {set charge [$sel get charge]}
foreach null {0} {set pos [$sel get z]}
set out [open $outputFile w]

# Get the force in piconewtons for this frame.
for {set i 0} {$i < $samplesZ} {incr i} {
    set rz [expr $zStart + ($zEnd-$zStart)/$samplesZ*$i]
    set fz 0.0

    foreach z $pos q $charge {
	set zCurr [expr $z + $rz]
	set ind [expr int(floor(($zCurr-$z0)/$dz))]
	if {[expr $ind >= 0 && $ind < $nz-1]} {
	    
	    set ezA [lindex $eFieldZ $ind]
	    set ezB [lindex $eFieldZ [expr $ind+1]]
	    set zA [expr $z0 + $ind*$dz]
	    set ez [expr $ezA + ($ezB-$ezA)/$dz*($zCurr-$zA)]
	    set fz [expr $fz + $q*$ez]
	}
    }
    set fz [expr 1602.176487*$fz]
    # Write the force in pN and the position in nm.
    puts "electrostatic force: [expr 0.1*$rz] $fz"
    puts $out "[expr 0.1*$rz] $fz"
}
close $out
$sel delete
exit



