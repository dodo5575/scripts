# Randomly distribute lipids in a box of a given size.
# to use: vmd -dispdev text -e randomlipbox.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>
package require psfgen

# Keep "nLipids" under 46656 if you know what's good for you.
set nLipids 400
set distance 3.0
set attempts 100

# Input:
set lipidpdb lipidcg.pdb
#Output:
set finalpdb lipidbox.pdb
set finalpsf lipidbox.psf

#Specify the box size.
set bx 100
set by 100
set bz 100
# Read the topology.
resetpsf
topology di21.top
topology lipid2.top

# Create a system with "nLipids" replicas of the lipid.
set base36 "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ"
for {set i 0} {$i < $nLipids} {incr i} {
	set c0 [string index $base36 [expr $i%36]]
	set c1 [string index $base36 [expr ($i/36)%36]]
	set c2 [string index $base36 [expr ($i/36/36)%36]]
	set seg [format "L%s%s%s" $c2 $c1 $c0]
	
	# Create a new segment for each lipid.
	segment $seg {
		first NONE
		last NONE
		pdb $lipidpdb
	}
	coordpdb $lipidpdb $seg
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
set unique_seg [lsort -unique [$all get segid]]
foreach seg $unique_seg {
    set sel [atomselect top [format "segid %s" $seg]]
    set pos0 [$sel get {x y z}]
    
    # Check that the lipids do not come too close.
    set violators 10 
    for {set i 0} {$violators > 0 && $i < $attempts} {incr i} {
	# Return the atoms to their original position.
	$sel set {x y z} $pos0
	
	# Calculate random angles and positions.
	set xt [expr (rand()-0.5)*$bx]
	set yt [expr (rand()-0.5)*$by]
	set zt [expr (rand()-0.5)*$bz]
	set xa [expr rand()*180]
	set ya [expr rand()*180]
	set za [expr rand()*180]

	# Put the lipid in a random orientation and position.
	$sel move [transaxis x $xa]
	$sel move [transaxis y $ya]
	$sel move [transaxis z $za]
    	$sel moveby "$xt $yt $zt"
	
	set vio [atomselect top \
		[format "segname %s and within %s of not segname %s" \
			$seg $distance $seg]]
	set violators [$vio num]
    }
    
    if {$i == $attempts} {
    	puts "WARNING: Too tired to fix lipid-lipid collision"
    }
}

# Write the final system.
$all writepdb $finalpdb
exit



