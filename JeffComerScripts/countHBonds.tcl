proc countHBonds {resid cutoff} {
	set oxygen "(resid $resid and name OH2)"
	set hydrogen1 "(resid $resid and name H1)"
	set hydrogen2 "(resid $resid and name H2)"
	
	set otherHydrogen "(name H1 H2 and not resid $resid)"
	set otherOxygen "(name OH2 and not resid $resid)"

	set o [atomselect top $oxygen]
	set oId [$o get {name resid}]
	set h1 [atomselect top $hydrogen1]
	set h1Id [$h1 get {name resid}]
	set h2 [atomselect top $hydrogen2]
	set h2Id [$h2 get {name resid}]
	
	set bo [atomselect top "$otherHydrogen and within $cutoff of $oxygen"]
	set bh1 [atomselect top "$otherOxygen and within $cutoff of $hydrogen1"]
	set bh2 [atomselect top "$otherOxygen and within $cutoff of $hydrogen2"]

	puts "Atom $oId is bound to:"
	puts [$bo get {name resid}]
	puts ""
	puts "Atom $h1Id is bound to:"
	puts [$bh1 get {name resid}]
	puts ""
	puts "Atom $h2Id is bound to:"
	puts [$bh2 get {name resid}]
	puts ""
	
	$o delete
	$h1 delete
	$h2 delete
	$bo delete
	$bh1 delete
	$bh2 delete
}



