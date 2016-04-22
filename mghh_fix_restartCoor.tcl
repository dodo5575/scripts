# This script fix the broken structure of mghh.
# Usage: vmd -dispdev text -e $argv0 -args inPSF inCoor outCoor  
# Author: Chen-Yu Li <cli56@illinois.edu>
# 2014/2/27

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

proc mghhFix {inPSF inCoor outCoor} {

    set OHA_t [list -1.794  0.738 -0.161] 
    set H1A_t [list -2.558  0.352  0.267] 
    set H2A_t [list -2.158  1.379 -0.771] 
    set OHB_t [list  0.064 -0.688 -1.861] 
    set H1B_t [list  0.544 -1.476 -2.115] 
    set H2B_t [list -0.442 -0.453 -2.639] 
    set OHC_t [list -0.773 -1.634  0.695] 
    set H1C_t [list -0.788 -1.968  1.592] 
    set H2C_t [list -1.187 -2.323  0.176] 
    set OHD_t [list  0.648  1.657 -0.733] 
    set H1D_t [list  0.220  2.495 -0.911] 
    set H2D_t [list  1.582  1.839 -0.830] 
    set OHE_t [list  0.102  0.755  1.712] 
    set H1E_t [list  0.625  1.493  2.025] 
    set H2E_t [list -0.607  0.677  2.350] 
    set OHF_t [list  1.826 -0.672  0.183] 
    set H1F_t [list  2.541 -0.208  0.618] 
    set H2F_t [list  1.985 -1.596  0.379] 
    
    set defaultD 2.3313897300712814
    
    set molID [mol load psf $inPSF namdbin $inCoor]
    
    set OHA [atomselect top "resname MGH and name OHA"]
    
    set resIDlist [$OHA get resid]
    
    foreach id $resIDlist {
    
        set MG [atomselect top "resid $id and name MG"]
    
        set MG_p [measure center $MG]
    
        set OHA [atomselect top "resid $id and name OHA"]
        $OHA moveby [vecSub [vecAdd $OHA_t $MG_p] [measure center $OHA]] 
        set H1A [atomselect top "resid $id and name H1A"]
        $H1A moveby [vecSub [vecAdd $H1A_t $MG_p] [measure center $H1A]] 
        set H2A [atomselect top "resid $id and name H2A"]
        $H2A moveby [vecSub [vecAdd $H2A_t $MG_p] [measure center $H2A]] 
    
        set OHB [atomselect top "resid $id and name OHB"]
        $OHB moveby [vecSub [vecAdd $OHB_t $MG_p] [measure center $OHB]] 
        set H1B [atomselect top "resid $id and name H1B"]
        $H1B moveby [vecSub [vecAdd $H1B_t $MG_p] [measure center $H1B]] 
        set H2B [atomselect top "resid $id and name H2B"]
        $H2B moveby [vecSub [vecAdd $H2B_t $MG_p] [measure center $H2B]] 
    
        set OHC [atomselect top "resid $id and name OHC"]
        $OHC moveby [vecSub [vecAdd $OHC_t $MG_p] [measure center $OHC]] 
        set H1C [atomselect top "resid $id and name H1C"]
        $H1C moveby [vecSub [vecAdd $H1C_t $MG_p] [measure center $H1C]] 
        set H2C [atomselect top "resid $id and name H2C"]
        $H2C moveby [vecSub [vecAdd $H2C_t $MG_p] [measure center $H2C]] 
    
        set OHD [atomselect top "resid $id and name OHD"]
        $OHD moveby [vecSub [vecAdd $OHD_t $MG_p] [measure center $OHD]] 
        set H1D [atomselect top "resid $id and name H1D"]
        $H1D moveby [vecSub [vecAdd $H1D_t $MG_p] [measure center $H1D]] 
        set H2D [atomselect top "resid $id and name H2D"]
        $H2D moveby [vecSub [vecAdd $H2D_t $MG_p] [measure center $H2D]] 
    
        set OHE [atomselect top "resid $id and name OHE"]
        $OHE moveby [vecSub [vecAdd $OHE_t $MG_p] [measure center $OHE]] 
        set H1E [atomselect top "resid $id and name H1E"]
        $H1E moveby [vecSub [vecAdd $H1E_t $MG_p] [measure center $H1E]] 
        set H2E [atomselect top "resid $id and name H2E"]
        $H2E moveby [vecSub [vecAdd $H2E_t $MG_p] [measure center $H2E]] 
    
        set OHF [atomselect top "resid $id and name OHF"]
        $OHF moveby [vecSub [vecAdd $OHF_t $MG_p] [measure center $OHF]] 
        set H1F [atomselect top "resid $id and name H1F"]
        $H1F moveby [vecSub [vecAdd $H1F_t $MG_p] [measure center $H1F]] 
        set H2F [atomselect top "resid $id and name H2F"]
        $H2F moveby [vecSub [vecAdd $H2F_t $MG_p] [measure center $H2F]] 
    
    #    foreach sym {A B C D E F} { 
    #        set water [atomselect top "resid $id and name OH${sym} H1${sym} H2${sym}"]
    #        set displacement [vecSub [measure center $MG] [measure center $water]]
    #        set distance [vecLength $displacement]
    #        if {$distance > 5} {
    #	    puts "ID: $id ; $sym ; $distance"
    #
    #	    $water moveby [vecScale [expr $distance - $defaultD] [vecScale [expr 1 / ($distance)] $displacement]]
    #        }
    #    }
    
    }

    set all [atomselect top all]
    animate write namdbin $outCoor
}

if {$argc < 2} {
    puts "vmd -dispdev text -e $argv0 -args inPSF inCoor outCoor"
    exit
}

set inPSF [lindex $argv 0]
set inCoor [lindex $argv 1]
set outCoor [lindex $argv 2]

mghhFix $inPSF $inCoor $outCoor 

exit
