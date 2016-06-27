
for {set i 0} {$i < 20} {incr i} {

    mol new ../triangle.psf
    mol addfile ../triangle.pdb

    foreach seg "T0 T1 T2 T3" {

        set sel [atomselect top "segname $seg"]
        $sel set segname ${seg}${i}
        $sel delete
    }


    set all [atomselect top all]


    $all move [transaxis x [expr rand() * 360]]
    $all move [transaxis y [expr rand() * 360]]
    $all move [transaxis z [expr rand() * 360]]

    $all moveby " [expr rand() * 500] 0 0 " 
    $all moveby " 0 [expr rand() * 500] 0 "
    $all moveby " 0 0 [expr rand() * 500] "

    $all writepsf triangle_${i}.psf
    $all writepdb triangle_${i}.pdb

    $all delete
    mol delete all

}


