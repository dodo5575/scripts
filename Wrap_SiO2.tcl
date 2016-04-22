proc Wrap {psfPrefix dcdPrefix outPut} {

    package require pbctools

    # Input:
    set psf $psfPrefix.psf
    set dcd $dcdPrefix.dcd

    # load trajectory and wait until loading completes 
    mol load psf $psf dcd $dcd
    set n [ molinfo top get numframes ]
    puts "  Loaded $n frames."

#    #pbc unwrap -all
#    
#    # set selections
#    set sel [atomselect top $selText]
#    set all [atomselect top all]
#
#    # define first and last frames and stride
#    set firstframe 0
#    set lastframe $n
#    set stride 1
#
#    $all frame 0
#    $sel frame 0    
#
#    # proceed frame by frame
#    for {set i $firstframe} {$i <= $lastframe} {set i [ expr $i + $stride ] } {
#    
#        $all frame $i
#        $sel frame $i    
#	 
#        # align dynamic selection to static one
#        set alignVector [measure center $sel]
#    
#        $all moveby [vecinvert $alignVector]
#        
#        puts "frame $i" 
#    }

    pbc wrap -all -sel "not (nucleic or segname U0)" -compound residue -center origin
    

    animate write dcd ${outPut}.dcd waitfor all top

    mol delete top

}


if {$argc < 1} {
    puts "vmd -dispdev text -e $argv0 -args psfPrefix dcdPrefix outPutPrefix"
    exit
}

set psfPrefix [lindex $argv 0]
set dcdPrefix [lindex $argv 1]
set outPut [lindex $argv 2]

Wrap $psfPrefix $dcdPrefix $outPut

exit

