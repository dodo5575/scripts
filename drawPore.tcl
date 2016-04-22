# This script will remove water from psf and pdf files.
# Use with: vmd -e drawPore.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

# Parameters:
proc drawPore {molid radius0 radius1 length {res 20}} {
#        graphics $molid delete all
#        graphics $molid color silver
	set segments $res
	set dz [expr double($length)/$res]
	
        set nz [expr 0.5*int(ceil($length/$dz))]
        # set zi [expr -0.5*$length]
        set zi 0
	set pi [expr 4.0*atan(1.0)]
	set phi [expr 2.0*$pi/$segments]

	for {set j 0} {$j < $nz} {incr j} {
		set z0 [expr $zi + $dz*$j]
		set z1 [expr $zi + $dz*($j+1.0)]
		set s0 [expr $radius0 + 2.0*($radius1-$radius0)/$length*abs($z0)]
		set s1 [expr $radius0 + 2.0*($radius1-$radius0)/$length*abs($z1)]
	
		for {set i 0} {$i < $segments} {incr i} {
			set rUnit0 [list [expr cos(($i-1)*$phi)] [expr sin(($i-1)*$phi)] 0.0]
			set rUnit1 [list [expr cos($i*$phi)] [expr sin($i*$phi)] 0.0]
		
			set ra [vecadd [vecscale $s0 $rUnit0] [list 0.0 0.0 $z0]]
			set rb [vecadd [vecscale $s1 $rUnit0] [list 0.0 0.0 $z1]]
			set rc [vecadd [vecscale $s1 $rUnit1] [list 0.0 0.0 $z1]]
			set rd [vecadd [vecscale $s0 $rUnit1] [list 0.0 0.0 $z0]]
		
			graphics $molid triangle $ra $rb $rc
			graphics $molid triangle $ra $rc $rd
		}
    # display update 
	}
}






