# computeForceOnProtein
#
# NOTE: this assumes a harmonic potential, i.e. consexp = 2
#
# author: dbwells2@uiuc.edu

if { $argc != 5 } {
    puts "Usage: vmd -dispdev text -e computeForceOnProtein.tcl -args <psf> <pdb> <dcd> <k>"
    exit
}

set psf	[lindex $argv 0]
set pdb	[lindex $argv 1]
set dcd	[lindex $argv 2]
set k	[expr 2 * [lindex $argv 3]]	;# because NAMD calculates V = kr^2 (or whatever consexp is) => F = 2kr

set refmol [mol load psf $psf pdb $pdb]
set tramol [mol load psf $psf dcd $dcd]

set refsel [atomselect $refmol "beta > 0"]
set trasel [atomselect $tramol "index [$refsel get index]"]

set refpos [$refsel get {x y z}]

set nframes [molinfo top get numframes]

set outCh [open ${dcd}.proteinForce "w"]

for { set frame 0 } { $frame < $nframes } { incr frame } {
    $trasel frame $frame
    set trapos [$trasel get {x y z}]
    
    set totalForce {0 0 0}
    foreach tpos $trapos rpos $refpos {
	set dpos [vecsub $rpos $tpos]
	set force_kcal [vecscale $dpos $k]
	set force [vecscale $force_kcal 69.523]
	
	set totalForce [vecadd $totalForce $force]
    }
    
    puts $outCh "$frame $totalForce"
    puts [format "FRAME %5s/%-5s  force = %s pN" $frame [expr { $nframes - 1 }] $totalForce]
}

close $outCh

exit
