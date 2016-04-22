# Auth: RCarr 10/8/2008
#
# These are common routines.
#

set kB_JperK 1.3806505e-23
set JtoKCalperMol 1.4393264e20
set dalton_per_kilogram 6.02214179e26
set PI 3.14159265358979323846

proc totalframes { dcdFile } {
    return [exec /home/rccarr2/bin/totalframes $dcdFile]
}

proc catdcd { outFile args } {
#    puts "This is the output $outFile"

#    puts "These are the args: $args"

    puts [lindex $args 0]
    eval "exec" "catdcd" "-o" $outFile [lindex $args 0]
    puts "Done with eval"
}

proc check_system { psf pdb } {
    # Calculate System Size
    set molID [mol load psf $psf pdb $pdb]
    
    set all [atomselect $molID "all"]
    set list [measure minmax $all]
    
    set size [vecsub [lindex $list 1] [lindex $list 0]]
    puts "\nSystem Size is $size"
    
    set q [measure sumweights $all weight charge]
    puts "Total charge: $q"
    
    mol delete $molID
}

# This is for 295 K
proc get_concentration { solute_seltext { molID -1 } } {

    if { $molID == -1 } {
	set molID [molinfo top]
    }

    set solute [atomselect $molID $solute_seltext]
    set water [atomselect $molID "water and noh"]

    set concentration [expr 55.455955 * double([$solute num]) / double ([$water num])]

    return $concentration
    
}

# use lsort - command random_sort
proc random_sort { A B } {
    if { [expr rand()] < 0.5 } {
	return -1
    } else {
	return 1
    }
}

# Aligns the v1 axis of a molecule along an axis
proc align_to_axis {sel v1 {axis z} {random 0}} {
    set v1 [vecnorm $v1]

    $sel move [transvecinv $v1]
    $sel moveby [vecinvert [measure center $sel]]

    if { $random } {
	$sel move [transabout [list 1 0 0] [expr {rand()*360}] deg]
	$sel moveby [vecinvert [measure center $sel]]
    }
    
    if { [string equal $axis z] } {
	set matrix [transaxis y -90 deg]
    } elseif { [string equal $axis -z] } {
	set matrix [transaxis y 90 deg]
    } elseif { [string equal $axis x] } {
	set matrix [transidentity]
    } elseif { [string equal $axis -x] } {
	set matrix [transaxis y [expr 180] deg]
    } elseif { [string equal $axis y] } {
	set matrix [transaxis z 90 deg]
    } elseif { [string equal $axis -y] } {
	set matrix [transaxis z -90 deg]
    } else {
	print "Error: align_to_axis accepts axes of (x,-x,y,-y,z,-z)"
    }

    $sel move $matrix

}


# Give this a vector containing the x y and z length of the system
# compatible with lindex \[pbc get \] 0
# Currently suppors only cubic systems.
proc random_position { system_vectors {center 1} } {

    if { [lindex $system_vectors 0] == 0 && [lindex $system_vectors 1] == 0 && [lindex $system_vectors 2] == 0 } {
	print "random_position: system vector is 0.0 0.0 0.0 !!!"
	return -1
    }


    set xyz [veczero]

    if { $center } {
	set xyz [list [expr { (rand()-0.5) * [lindex $system_vectors 0]}] [expr {(rand()-0.5) * [lindex $system_vectors 1]}]  [expr {(rand()-0.5) * [lindex $system_vectors 2]}] ]
    } else {
	set xyz [list [expr {rand() * [lindex $system_vectors 0]}] [expr {rand() * [lindex $system_vectors 1]}]  [expr {rand() * [lindex $system_vectors 2]}] ]
    }

    return $xyz
}


proc recolor_dx { {molID -1} } {

    if { $molID == -1 } {
	set molID [molinfo top]
    }
    
    set minmax [mol scaleminmax $molID 0]
    color scale min 0.0
    color scale max 1.0
    color scale midpoint [expr {abs([lindex $minmax 0]) / ([lindex $minmax 1]-[lindex $minmax 0])}]
    
}

proc load_dx { dxFile } {
    set molID [mol new $dxFile]
    mol delrep 0 $molID

    # The bottom slab
    mol color Volume 0 
    mol representation VolumeSlice 0.0 0.0 0.0 2.0
    mol addrep $molID
    mol representation VolumeSlice 0.0 0.0 1.0 2.0
    mol addrep $molID

    recolor_dx $molID

    rotate x by -90
    
    return $molID
}