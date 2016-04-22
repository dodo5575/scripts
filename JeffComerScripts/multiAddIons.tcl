# Add ions (KCl) of given ionic strength.
# vmd -dispdev text -e addIons.tcl

set concList {100mM 300mM 1000mM 1300mM 1700mM}
set numList {7 21 69 90 118}
set alreadyNum 5

proc addIons {conc number} {
    # Parameters:
    # Define which ions to replace with which ions.
    set ionFrom "SOD"
    set ionTo "POT"
    set name pore_none_basepair
    # Input:
    set topology top_all27_prot_lipid_pot.inp
    set psf ${name}1.psf
    set pdb ${name}1.pdb
    # Output:
    set prefixIons ${name}_NaCl
    set prefix bulk${conc}

    package require autoionize
    autoionize -psf $psf -pdb $pdb -ncl $number -nna $number -o $prefixIons -from 2.0 -between 2.0

    set psfFile ${prefixIons}.psf
    set pdbFile ${prefixIons}.pdb

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
    mol delete top
}

foreach conc $concList number $numList {
    addIons $conc [expr {$number-$alreadyNum}] 
}