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

proc mghhFix {inPrefix outPrefix} {

    set molID [mol load psf ${inPrefix}.psf pdb ${inPrefix}.pdb]
    
    set H2O [atomselect top water]

    foreach {x0 y0 z0} [lindex [measure minmax $H2O] 0] {}
    foreach {x1 y1 z1} [lindex [measure minmax $H2O] 1] {}
    set xL [expr $x1 - $x0] 
    set yL [expr $y1 - $y0] 
    set zL [expr $z1 - $z0] 
   
    set mg1 [atomselect top "resname MGH and name MG"]

    set segname [lsort -unique [$mg1 get segname]] 
    foreach seg $segname {

        set mg2 [atomselect top "segname $seg and name MG"]
        set resID [$mg2 get resid] 
 
        set fix_list ""
        foreach id $resID {

            set mg3 [atomselect top "segname $seg and resid $id and name MG"]
            set mgCoor [measure center $mg3]        

            set mgh [atomselect top "segname $seg and resid $id and (not name MG)"]
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

                $sel delete
            }
            $mgh delete
            $mg3 delete
        }
        $mg2 delete

    }
    set all [atomselect top all]
    $all writepsf ${outPrefix}.psf
    $all writepdb ${outPrefix}.pdb
}

if {$argc < 2} {
    puts "vmd -dispdev text -e $argv0 -args inPrefix outPrefix"
    exit
}

set inPrefix  [lindex $argv 0]
set outPrefix [lindex $argv 1]

mghhFix $inPrefix $outPrefix 

exit
