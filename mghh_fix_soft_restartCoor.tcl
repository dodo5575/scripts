# This script fix the broken structure of mghh.
# Usage: vmd -dispdev text -e $argv0 -args inPrefix outPrefix  
# Author: Chen-Yu Li <cli56@illinois.edu>
# 2015/7/8

proc vecAdd {a b} {
    set c {}
    foreach ai $a bi $b {
        lappend c [expr $ai+$bi]
    }
    return $c
}

proc vecSub {a b} {
    set c {}
    foreach ai $a bi $b {
        lappend c [expr $ai-$bi]
    }
    return $c
}

proc vecLength {a} {
    set sum 0
    foreach ai $a {
        set sum [expr $sum + $ai*$ai]
    }
    return [expr sqrt($sum)]
}

proc vecScale {s a} {
    set b {}
    foreach ai $a {
        lappend b [expr $ai*$s]
    }
    return $b
}

proc mghhFix {inPsf inCoor outCoor} {

    #set OHA_t [list -1.794  0.738 -0.161] 
    #set H1A_t [list -2.558  0.352  0.267] 
    #set H2A_t [list -2.158  1.379 -0.771] 
    #set OHB_t [list  0.064 -0.688 -1.861] 
    #set H1B_t [list  0.544 -1.476 -2.115] 
    #set H2B_t [list -0.442 -0.453 -2.639] 
    #set OHC_t [list -0.773 -1.634  0.695] 
    #set H1C_t [list -0.788 -1.968  1.592] 
    #set H2C_t [list -1.187 -2.323  0.176] 
    #set OHD_t [list  0.648  1.657 -0.733] 
    #set H1D_t [list  0.220  2.495 -0.911] 
    #set H2D_t [list  1.582  1.839 -0.830] 
    #set OHE_t [list  0.102  0.755  1.712] 
    #set H1E_t [list  0.625  1.493  2.025] 
    #set H2E_t [list -0.607  0.677  2.350] 
    #set OHF_t [list  1.826 -0.672  0.183] 
    #set H1F_t [list  2.541 -0.208  0.618] 
    #set H2F_t [list  1.985 -1.596  0.379] 
    #
    #set defaultD 2.3313897300712814
    #
    set molID [mol load psf $inPsf namdbin $inCoor]
    
    set H2O [atomselect top water]

    foreach {x0 y0 z0} [lindex [measure minmax $H2O] 0] {}
    foreach {x1 y1 z1} [lindex [measure minmax $H2O] 1] {}
    set xL [expr $x1 - $x0] 
    set yL [expr $y1 - $y0] 
    set zL [expr $z1 - $z0] 
   
    set mg [atomselect top "resname MGH and name MG"]

    set resID [$mg get resid] 
 
    set fix_list ""
    foreach id $resID {

        set mg [atomselect top "resname MGH and resid $id and name MG"]
        set mgCoor [measure center $mg]        

        set mgh [atomselect top "resname MGH and resid $id and (not name MG)"]
        set indices [$mgh get index]

        foreach ind $indices {
            set sel [atomselect top "index $ind"]
            set coor [measure center $sel]
            if {[expr [lindex $coor 0] - [lindex $mgCoor 0] ] > 5} {$sel moveby "[expr -1 * $xL] 0 0"} 
            if {[expr [lindex $coor 0] - [lindex $mgCoor 0] ] <-5} {$sel moveby "[expr  1 * $xL] 0 0"} 
            if {[expr [lindex $coor 1] - [lindex $mgCoor 1] ] > 5} {$sel moveby "0 [expr -1 * $yL] 0"}
            if {[expr [lindex $coor 1] - [lindex $mgCoor 1] ] <-5} {$sel moveby "0 [expr  1 * $yL] 0"}
            if {[expr [lindex $coor 2] - [lindex $mgCoor 2] ] > 5} {$sel moveby "0 0 [expr -1 * $zL]"}
            if {[expr [lindex $coor 2] - [lindex $mgCoor 2] ] <-5} {$sel moveby "0 0 [expr  1 * $zL]"}
        }

    }

    set all [atomselect top all]
    animate write namdbin $outCoor
}

if {$argc < 2} {
    puts "vmd -dispdev text -e $argv0 -args inPsf inCoor outCoor"
    exit
}

set inPsf   [lindex $argv 0]
set inCoor  [lindex $argv 1]
set outCoor [lindex $argv 2]

mghhFix $inPsf $inCoor $outCoor 

exit
