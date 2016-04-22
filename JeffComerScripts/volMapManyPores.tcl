# This script will remove all residues in the selection
# from psf and pdf files.
# Use with: vmd -dispdev text -e removeResidues.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set radMax 20
set zMax 20
set selTextGeo "abs(z) < $zMax and abs(x) < $radMax and abs(y) < $radMax"
set selTextSi "resname SIN and name \"SI.*\""
set selTextN "resname SIN and name \"N.*\""
set halfRadMinSi [expr 2.135]
set halfRadMinN [expr 1.9975]
set res 0.4
# Input:
set nameList {pore1.2 pore1.6 pore2.0 pore2.0_thick pore2.0_trap pore2.2}
set psfList {pore1.2_all.psf pore1.6_coil_expedite.psf pore2.0_all.psf pore+dna_E4V.psf trap2.0.psf pore2.2_trap.psf}
set coorList {coil_p1.2_4Vb1.restart.coor coilx_p1.6_4V14.restart.coor coil_p2.0_2V4.restart.coor coil_cont_1V8.restart.coor trap2.0_0V3.restart.coor pore2.2_trap_eq2.restart.coor}
# Output:
set outPrefix dist_big_

for {set pore 0} {$pore < [llength $nameList]} {incr pore} {
    set name [lindex $nameList $pore]
    set psf [lindex $psfList $pore]
    set coor [lindex $coorList $pore]

    mol load psf $psf
    mol addfile $coor
    set sel [atomselect top $selTextSi]
    $sel set radius $halfRadMinSi
    $sel delete
    set sel [atomselect top $selTextN]
    $sel set radius $halfRadMinN
    $sel delete

    set sel [atomselect top "(($selTextSi) or ($selTextN)) and ($selTextGeo)"]
    puts "Generating volmap from [$sel num] atoms."
    volmap distance $sel -o ${outPrefix}_${name}.dx -res $res -cutoff 5.0
    $sel delete

    mol delete top
}
exit



