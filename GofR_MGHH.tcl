# This script calculate the g(r) in MGHH.
# Usage: vmd -dispdev text -e $argv0 -args name PSF outDir COOR 
# Modified by Chen-Yu Li <cli56@illinois.edu>
# 2014.5.25


proc vecLength {a} {
    set sum 0
    foreach ai $a {
        set sum [expr $sum + $ai*$ai]
    }
    return [expr sqrt($sum)]
}
proc vecSub {a b} {
    set c {}
    foreach ai $a bi $b {
        lappend c [expr $ai-$bi]
    }
    return $c
}


proc GofR {name PSF outDir COOR} {

    set min 0
    set max 10
    set binSize 0.1
    set binN [expr int(($max - $min) / $binSize)]
    #puts $binN

    array set distListIndex {}
    array set distList {}
    for {set i 0} {$i < $binN} {incr i} {
        set distListIndex($i) [expr $min + (0.5 + $i) * $binSize]
        set distList($i) 0
    }

    for {set i 0} {$i < [array size distListIndex]} {incr i} {
        puts $distListIndex($i)
    }
    for {set i 0} {$i < [array size distList]} {incr i} {
        puts $distList($i)
    }
    #puts $distListIndex
    #puts $distList
    #puts [llength $distListIndex]

    mol load psf $PSF namdbin $COOR

    set sel1 [atomselect top "segname MGHH and name MG"]
    set sel2 [atomselect top "segname MGHH and name OHA OHB OHC OHD OHE OHF"]

    set index1 [$sel1 get index]    
    set sel1N [llength $index1]
    set index2 [$sel2 get index]    

    foreach ind1 $index1 {

        set pos1 [measure center [atomselect top "index $ind1"]]

        foreach ind2 $index2 {

            set pos2 [measure center [atomselect top "index $ind2"]]

            set dist [vecLength [vecSub $pos1 $pos2]]
            
            if {$dist >= $max} {
                continue
            } else {

                set d [expr int($dist / $binSize)]
                set distList($d) [expr $distList($d) + 1]

            }

        }

    }

    #for {set i 0} {$i < [array size distList]} {incr i} {
    #    puts $distList($i)
    #}

    set out [open "${outDir}/${name}.dat" w]
    for {set i 0} {$i < [array size distList]} {incr i} {

        #set distList [lreplace $distList $i $i [expr [lindex $distList $i] / $sel1N ]]
        set interval $distListIndex($i)
        set dist  [expr $distList($i)/$sel1N]
        puts $out "$interval\t$dist"

    }


}


if {$argc < 4} {
    puts "vmd -dispdev text -e $argv0 -args name PSF outDir COOR"
    exit
}

set name [lindex $argv 0]
set PSF [lindex $argv 1]
set outDir [lindex $argv 2]
set COOR [lindex $argv 3]

GofR $name $PSF $outDir $COOR

exit

