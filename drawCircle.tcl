# This script draw circles at given position.
# Use with: vmd -e drawCircle.tcl
# Author: Chen-Yu Li <cli56@illinois.edu>
# 2014/10/27


# Parameters:
proc drawCircle {molid x y z radius0 radius1 {res 20}} {
#        graphics $molid delete all
#        graphics $molid color silver
	set segments $res
	#set dz [expr double($length)/$res]
	
        #set nz [expr 0.5*int(ceil($length/$dz))]
        set dr [expr ($radius1 - $radius0)/$segments]
        set nr [expr int(ceil(($radius1 - $radius0)/$dr))]
        # set zi [expr -0.5*$length]
        set zi 0
	set pi [expr 4.0*atan(1.0)]
	set phi [expr 2.0*$pi/$segments]

	for {set j 0} {$j < $nr} {incr j} {
		#set z0 [expr $zi + $dz*$j]
		#set z1 [expr $zi + $dz*($j+1.0)]
		set s0 [expr $radius0 + $dr*$j]
		set s1 [expr $radius0 + $dr*($j+1)]
	
		for {set i 0} {$i < $segments} {incr i} {
			#set rUnit0 [list [expr cos(($i-1)*$phi) + $x] [expr sin(($i-1)*$phi) + $y] $z]
			#set rUnit1 [list [expr cos($i*$phi) + $x] [expr sin($i*$phi) + $y] $z]
			set rUnit0 [list [expr cos(($i-1)*$phi)] [expr sin(($i-1)*$phi)] 0.0]
			set rUnit1 [list [expr cos($i*$phi)] [expr sin($i*$phi)] 0.0]
		
			set ra [vecadd [list $x $y $z] [vecscale $s0 $rUnit0]] 
			set rb [vecadd [list $x $y $z] [vecscale $s1 $rUnit0]] 
			set rc [vecadd [list $x $y $z] [vecscale $s1 $rUnit1]] 
			set rd [vecadd [list $x $y $z] [vecscale $s0 $rUnit1]]
		
			graphics $molid triangle $ra $rb $rc
			graphics $molid triangle $ra $rc $rd
		}
    # display update 
	}
}






