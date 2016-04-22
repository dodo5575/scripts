proc quiet {} {}

source /home/jcomer/scripts/vector.tcl

mol load psf mspa_ions.psf pdb mspa_ions.pdb
set sel [atomselect top "ions"]
set posList [$sel get {x y z}]; quiet

mol load psf urea_good_tiled.psf pdb urea_good_tiled.pdb
set all [atomselect top "all"]
set segList [$all get segname]; quiet

set ind 0
foreach pos $posList {
    set seg [lindex $segList $ind]

    set s [atomselect top "segname $seg"] 
    set cen [measure center $s weight mass]
    $s moveby [vecinvert $cen]
    $s move [matMake4 [matRandomRot]]
    $s moveby $pos
    $s delete

    incr ind
}
puts "Did $ind"

for {set i $ind} {$i < [llength $segList]} {incr i} {
    set seg [lindex $segList $i]
    set s [atomselect top "segname $seg"]
    set cen [measure center $s weight mass]
    $s moveby [vecinvert $cen]
    $s moveby {0 0 1000}
    $s delete
}

$all writepdb urea_scatter.pdb
$all delete
mol delete top
mol delete top

exit
