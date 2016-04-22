## usage: makePsf villin-gms.pdb  villin-noWater

proc makeDNAConstraints {inPrefix outPrefix} {
    set ID [mol new ${inPrefix}.psf]
    mol addfile ${inPrefix}.pdb $ID

    set all [atomselect $ID all]
    $all set occupancy 0.0
    $all set beta 0.0

    set atomsToConstrain [atomselect $ID "nucleic"]

    $atomsToConstrain set beta 1.0

    $all writepdb ${outPrefix}.pdb

    mol delete $ID

}

