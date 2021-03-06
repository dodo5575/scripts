# Add ions (KCl) of given ionic strength.
# vmd -dispdev text -e addIons.tcl

# Parameters:
set conc 0.1
# Define which ions to replace with which ions.
set ionFrom "SOD"
set ionTo "POT"
set name sys_cone
# Input:
set topology top_all27_prot_lipid_pot.inp
set psf ${name}_hex${sys}.psf
set pdb ${name}_hex${sys}.pdb
# Output:
set prefixIons ${name}_NaCl${sys}
set prefix bam_spec_pull${sys}

mol load psf $psf pdb $pdb

# Compute the number of ions to add.
set posSel [atomselect top "name POT"]
set posNum [$posSel num]
set negSel [atomselect top "name CLA"]
set negNum [$negSel num]
set sel [atomselect top "name OH2"]
set nw [$sel num]
set alpha 55.523
set ni [expr int(floor($conc*$nw/($alpha + 2.0*$conc)+0.5))]
$sel delete

# Get the charge.
set all [atomselect top all]
set other [atomselect top "not name POT CLA"]
set nq [expr int(floor([measure sumweights $other weight charge]+0.5))]
set nna [expr $ni - $posNum - $nq]
set ncl [expr $ni - $negNum]
puts "posNum0: $posNum"
puts "posNumQ: $nq"
puts "posNum1: $nna"
puts "negNum0: $negNum"
puts "negNum1: $ncl"
$all delete
$other delete
mol delete top

package require autoionize
#autoionize -psf $psf -pdb $pdb -is [expr 2.*$conc] -o $prefixIons -from 2.0 -between 2.0
autoionize -psf $psf -pdb $pdb -nna $nna -ncl $ncl -o $prefixIons -from 2.0 -between 2.0

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
exit

