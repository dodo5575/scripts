#Attention: use '_' instead of ' ' in selection text

proc AlignWrap {psfPrefix pdbPrefix dcdPrefix selText outPut} {

    package require pbctools

    # Input:
    set psf $psfPrefix.psf
    set pdb $psfPrefix.pdb
    set dcd $dcdPrefix.dcd

    mol load psf $psf pdb $pdb
    set sel [atomselect top $selText]
    set ref [measure center $sel]    


    # load trajectory and wait until loading completes 
    mol load psf $psf dcd $dcd
    set n [ molinfo top get numframes ]
    puts "  Loaded $n frames."

    #pbc unwrap -all
    
    # set selections
    set sel [atomselect top $selText]
    set all [atomselect top all]

    # define first and last frames and stride
    set firstframe 0
    set lastframe $n
    set stride 1

    $all frame 0
    $sel frame 0    

    # proceed frame by frame
    for {set i $firstframe} {$i <= $lastframe} {set i [ expr $i + $stride ] } {
    
        $all frame $i
        $sel frame $i    
	 
        # align dynamic selection to static one
        set alignVector [vecsub $ref [measure center $sel]]
    
        $all moveby $alignVector
        
        puts "frame $i" 
    }

    pbc wrap -all -sel "not ($selText or nucleic)" -compound residue -center origin
    

    animate write dcd ${outPut}.dcd waitfor all top

    mol delete top

}


if {$argc < 5} {
    puts "vmd -dispdev text -e $argv0 -args psfPrefix pdbPrefix dcdPrefix selText outPutPrefix"
    exit
}

set psfPrefix [lindex $argv 0]
set pdbPrefix [lindex $argv 1]
set dcdPrefix [lindex $argv 2]
#map the '_' in selection text to ' '
set selText   [lindex $argv 3]
set selText   [string map {_ \ } $selText]
set outPut    [lindex $argv 4]

AlignWrap $psfPrefix $pdbPrefix $dcdPrefix $selText $outPut

exit

