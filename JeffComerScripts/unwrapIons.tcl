set stride 1

package require pbctools

foreach sys {q-1 q2c} {
    mol load psf nw_anneal_dopc_${sys}.psf
    mol addfile nw_field0_dopc_a${sys}.dcd step $stride waitfor all
    
    set nFrames [molinfo top get numframes]
    puts [format "Reading %i frames." $nFrames]
    set last [expr {$nFrames-1}]
    
    pbc unwrap -all -sel "ions"
    puts "Unwrapped."

    set all [atomselect top all]
    animate write dcd nw_field0_dopc_${sys}_unwrap.dcd beg 0 end $last waitfor all sel $all top
    $all delete
    mol delete top
}
