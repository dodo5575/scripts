# Author: Jeff Comer <jcomer2@illinois.edu>

# Parameters:
set selText ions
set trailNum 16
set dirCol charge
set distDiffuse 0.1
set distDrift 0.3
# Input:
set psf nw_pore2.0_all.psf
set pdb nw_pore2.0_all.pdb
# Output:
set outName nw_traj

proc randomDirection {} {
    set pi [expr {4.0*atan(1.0)}]

    set cosPhi [expr {2.0*rand()-1.0}]
    set sinPhi [expr {sqrt(1.0-$cosPhi*$cosPhi)}]
    set theta [expr {2.0*rand()*$pi}]
    

    set x [expr {$sinPhi*cos($theta)}]
    set y [expr {$sinPhi*sin($theta)}]
    set z $cosPhi
    
    return [list $x $y $z]
}

mol load psf $psf pdb $pdb
set all [atomselect top all]
set sel [atomselect top $selText]
foreach quiet {0} {
    set charge [$sel get $dirCol]
    set pos [$sel get {x y z}]
}

set displace {}
foreach r $pos q $charge {
    set drift [list 0.0 0.0 [expr {$q*$distDrift}]]
    set diffuse [vecscale $distDiffuse [randomDirection]]
    lappend displace [vecadd $drift $diffuse]
}

for {set i 0} {$i < $trailNum} {incr i} {
    set newPos {}
    foreach r $pos d $displace {
	lappend newPos [vecadd $r [vecscale $i $d]]
    }
    
    $sel set {x y z} $newPos
    $all writepdb ${outName}${i}.pdb
}

mol delete top
$sel delete
$all delete
exit
