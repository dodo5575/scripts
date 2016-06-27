package require psfgen

set base   [lindex $argv 0]
set num    [lindex $argv 1]

resetpsf
for {set i 0} {$i < $num} {incr i} {

    set Prefix "${base}_${i}"

    readpsf  ${Prefix}.psf
    coordpdb ${Prefix}.pdb

}

writepsf [lindex $argv 2].psf
writepdb [lindex $argv 2].pdb

exit
