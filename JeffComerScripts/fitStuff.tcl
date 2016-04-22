mol load psf stp_CORE.psf pdb stp_CORE.pdb
mol load psf avidin_biotin_dna.psf pdb a-b_placement.pdb

set molList {0 1}
set segList {BTNA BTN}
set atomNames "SI O3I C10I N2I"

set selList {}
foreach m $molList seg $segList {
    lappend selList [atomselect $m "segname $segList and name $atomNames"]
}

set all [atomselect [lindex $molList 0] all]
#set pos [$all get {x y z}]
set mat [measure fit [lindex $selList 0] [lindex $selList 1] weight mass]

$all move $mat

$all writepdb strept_fit.pdb
exit
