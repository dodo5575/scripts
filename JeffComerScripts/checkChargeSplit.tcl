# Author: Jeff Comer <jcomer2@illinois.edu>

proc check {z} {
    set sel [atomselect top "z > $z"]
    set q [measure sumweights $sel weight charge]
    
    set typeList [lsort -unique [$sel get type]]
    foreach type $typeList {
	set s [atomselect top "type $type and z > $z"]
	puts "atom $type: [$s num]"
	$s delete
    } 

    $sel delete
    return $q
}
