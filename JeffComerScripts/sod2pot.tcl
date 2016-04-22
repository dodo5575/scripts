proc sod2pot {psfFile pdbFile topology prefix} {
    set ionFrom "SOD"
    set ionTo "POT"

    package require psfgen
    topology $topology

    puts "\nSod2pot) Reading ${psfFile}/${pdbFile}..."
    resetpsf
    readpsf $psfFile
    coordpdb $pdbFile
    mol load psf $psfFile pdb $pdbFile

    set sel [atomselect top "name $ionFrom"]
    foreach zero {0} {
	set poslist [$sel get {x y z}]
	set seglist [$sel get segid]
	set reslist [$sel get resid]
    }
    set num [llength $reslist]
    puts "Sod2pot) Found ${num} ${ionFrom} ions to replace..."

    set num 0
    foreach segid $seglist resid $reslist {
	delatom $segid $resid
	incr num
    }
    puts "Sod2pot) Deleted ${num} ${ionFrom} ions"

    segment $ionTo {
	first NONE
	last NONE
	foreach res $reslist {
	    residue $res $ionTo
	}
    }
    set num [llength $reslist]
    puts "Sod2pot) Created ${num} topology entries for ${ionTo} ions"

    set num 0
    foreach xyz $poslist res $reslist {
	coord $ionTo $res $ionTo $xyz
	incr num
    }
    puts "Sod2pot) Set coordinates for ${num} ${ionTo} ions"

    writepsf "${prefix}.psf"
    writepdb "${prefix}.pdb"
    puts "Sod2pot) Wrote ${prefix}.psf/${prefix}.pdb"
    puts "Sod2pot) All done."
    mol delete top
}
