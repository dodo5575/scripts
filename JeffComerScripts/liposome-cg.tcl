# Randomly distribute lipid pairs onto a sphere to form a liposome.
# to use: vmd -dispdev text -e randomlipbox.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>
package require psfgen

# Specify the liposome size.
set radius 80.0
# Keep "nLipids" under 46656 if you know what's good for you.
set nLipids 1800
set distance 1.0
set attempts 200

# Input:
# The lipid partner files should contain the lipid pair along the z-axis.
set lipid0pdb lipidpartner0.pdb
set lipid1pdb lipidpartner1.pdb
# Output:
set finalpdb lipo80.pdb
set finalpsf lipo80.psf
set temppdb tmp.pdb

# Read the topology.
topology lipid2.top

# Create a system with "nLipids" replicas of the lipid.
set base36 "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ"
for {set i 0} {$i < $nLipids} {incr i 2} {
	set c0 [string index $base36 [expr $i%36]]
	set c1 [string index $base36 [expr ($i/36)%36]]
	set c2 [string index $base36 [expr ($i/36/36)%36]]
	set seg [format "L%s%s%s" $c2 $c1 $c0]
	
	# Create a new segment for the first partner.
	segment $seg {
		first NONE
		last NONE
		pdb $lipid0pdb
	}
	coordpdb $lipid0pdb $seg
	
	set c0 [string index $base36 [expr ($i+1)%36]]
	set c1 [string index $base36 [expr (($i+1)/36)%36]]
	set c2 [string index $base36 [expr (($i+1)/36/36)%36]]
	set seg [format "L%s%s%s" $c2 $c1 $c0]
	
	# Create a new segment for the first partner.
	segment $seg {
		first NONE
		last NONE
		pdb $lipid1pdb
	}
	coordpdb $lipid1pdb $seg
	guesscoord
}

writepdb $finalpdb
writepsf $finalpsf

# Load the box of lipids.
mol load psf $finalpsf
mol addfile $finalpdb waitfor all
set all [atomselect top "all"]
set center [measure center $all weight mass]
$all moveby [vecinvert $center]

# Randomly distribute the lipids.
set tired 0
set count 0
set pi [expr 4.0*atan(1.0)]
set unique_seg [lsort -unique [$all get segid]]
foreach {seg0 seg1} $unique_seg {
    set sel [atomselect top [format "segid %s %s" $seg0 $seg1]]
    set pos0 [$sel get {x y z}]
    
    # Check that the lipids do not come too close.
    set violators 10 
    for {set i 0} {$violators > 0 && $i < $attempts} {incr i} {
	# Return the atoms to their original position.
	$sel set {x y z} $pos0
	
	# Create a random rotation matrix, uniform on a sphere.
	set psi [expr 2.0*rand()*$pi]
	set phi [expr acos(2.0*rand()-1.0)]
	set theta [expr 2.0*rand()*$pi]
	set basis [transaxis z $psi rad]
	set basis [transmult [transaxis x $phi rad] $basis]
	set basis [transmult [transaxis z $theta rad] $basis]
		
	# Rotate the lipid pairs and move to the surface of the sphere.
	$sel move $basis
	set z [vecnorm [coordtrans $basis {0.0 0.0 1.0}]]
	$sel moveby [vecscale $z $radius]
	
	if {!$tired} {
	set vio [atomselect top \
		[format "segname %s %s and within %s of not segname %s %s" \
		$seg0 $seg1 $distance $seg0 $seg1]]
	set violators [$vio num]
	} else { set violators 0 }
    }
    
    if {$i == $attempts} {
    	puts [format \
	"WARNING: Too tired to fix lipid-lipid collisions, %s lipids" \
	$count]
	set tired 1
    }
    incr count 2
}

# Write the final system.
$all writepdb $finalpdb
exit



