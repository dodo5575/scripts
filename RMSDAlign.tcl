

proc RMSDAlignWrap {psfPrefix pdbPrefix dcdPrefix selText outPut} {

    package require pbctools

    # Input:
    set psf $psfPrefix.psf
    set pdb $pdbPrefix.pdb
    set dcd $dcdPrefix.dcd

    # set reference selection
    mol load psf $psf pdb $pdb
    set selref [atomselect top $selText]
    puts "number of atoms in the selection: [$selref num]"

    # load trajectory and wait until loading completes 
    mol load psf $psf dcd $dcd
    set n [ molinfo top get numframes ]
    puts "  Loaded $n frames."

    #pbc unwrap -all
    
    # set dynamic selections
    set seldyn [atomselect top $selText]
    
    # define first and last frames and stride
    set firstframe 0
    set lastframe $n
    set stride 1
    #set selref [atomselect top $selection frame firstframe]
    
    # proceed frame by frame
    for {set i $firstframe} {$i <= $lastframe} {set i [ expr $i + $stride ] } {
    
        
        set seldyn [atomselect top $selText frame $i]
        set selAll [atomselect top all frame $i]
        # align dynamic selection to static one
        set matrix [ measure fit $seldyn $selref ]
    
        $selAll move $matrix
    
        # measure RMSD for all heavy atoms and for protein
        
        $selAll delete
        $seldyn delete
        
        puts "frame $i" 
    }

    #pbc wrap -all -sel "not nucleic" -compound residue -center origin
    

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
set selText [lindex $argv 3]
set outPut [lindex $argv 4]

RMSDAlignWrap $psfPrefix $pdbPrefix $dcdPrefix $selText $outPut

exit



