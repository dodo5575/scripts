namespace eval ::tclForces::forces {
    lappend exportProcs bond bondX bondAxis

    ## BOND: restrains the distance between 
    proc bond { ids1 ids2 k r0 {print 0} {energy_name ""} } {
	## calculate force
	foreach {c1 m1} [getCOM $ids1] { break }
	foreach {c2 m2} [getCOM $ids2] { break }
	set r [vecsub $c2 $c1]
	set rl [veclength $r] 
	set dr [expr {$rl-$r0}] ;# seperation distance

	## apply force
	#set force [vecscale $r [expr $dr*$k/([veclength $r])]]
	set force [vecscale $r [expr {$dr*$k/$rl}]]
	applyForce $ids1 $force masses $m1
	applyForce $ids2 $force masses -$m2
	
	## print force
	if { $print } { printForce bond seperation $dr force $force }
	
	if { [llength $energy_name] > 0 } { addEnergy $energy_name [expr {.5*$k*$dr*$dr}] }
    }
    proc bondX { ids1 ids2 k r0 {print 0} {energy_name ""} } {
	## calculate force
	foreach {c1 m1} [getCOM $ids1] { break }
	foreach {c2 m2} [getCOM $ids2] { break }
	set x1 [lindex $c1 0]
	set x2 [lindex $c2 0]
	set r [expr {$x2-$x1}]
	set dr [expr {abs($r)-$r0}] ;# seperation distance
	set sign [expr {$r/abs($r)}]

	## apply force
	set force "[expr {$sign*$dr*$k}] 0 0"
	applyForce $ids1 $force masses $m1
	applyForce $ids2 $force masses -$m2
	
	## print force
	if { $print } { printForce bondX dx $r ddx $dr force $force }
	if { [llength $energy_name] > 0 } { addEnergy $energy_name [expr {.5*$k*$dr*$dr}] }
    }
    
    proc bondAxis { ids1 ids2 k r0 axis {print 0} {energy_name ""} } {
	## calculate force
	foreach {c1 m1} [getCOM $ids1] { break }
	foreach {c2 m2} [getCOM $ids2] { break }
	set axis [vecnorm $axis]
	set r [vecproject [vecsub $c2 $c1] $axis]
	set rl [veclength $r]
	set sign [sign [vecdot $axis $r]]
	set dr [expr {$sign*$rl-$r0}] ;# seperation distance
	
	## apply force
	#set force [vecscale $r [expr $dr*$k/([veclength $r])]]
	set force [vecscale $axis [expr {$dr*$k}]]
	applyForce $ids1 $force masses $m1
	applyForce $ids2 $force masses -$m2
	
	## print force
	if { $print } { printForce bondAxis seperation $dr force $force }
	
	if { [llength $energy_name] > 0 } { addEnergy $energy_name [expr {.5*$k*$dr*$dr}] }
    }
}
