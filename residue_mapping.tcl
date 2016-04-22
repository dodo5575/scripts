# create a residue mapping between 2 pdbs
# Usage: vmd -dispdev text -e residue_mapping.tcl -args pdb1 pdb2 output 
# Author: Chen-Yu Li <cli56@illinois.edu> 
# 2015/7/23


set pdb1    [lindex $argv 0]
set pdb2    [lindex $argv 1]
set output  [lindex $argv 2]

set out [open $output w]
puts $out "$pdb1\t$pdb2"

set id1 [mol load pdb $pdb1]
set id2 [mol load pdb $pdb2]

set all1 [atomselect $id1 "all"]

set residue1 [lsort -unique [$all1 get residue]]

foreach r $residue1 {

    set sel1 [atomselect $id1 "residue $r"]

    set seg   [lsort -unique [$sel1 get segname]]
    set chain [lsort -unique [$sel1 get chain]]
    set resid [lsort -unique [$sel1 get resid]]

    set sel2 [atomselect $id2 "segname $seg and chain $chain and resid $resid"]

    set residue2 [lsort -unique [$sel2 get residue]]

    puts $out "$r\t$residue2"

    $sel1 delete
    $sel2 delete
}

close $out

