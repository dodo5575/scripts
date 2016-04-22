

proc RMSD {psfPrefix dcd lipidResname dim interval R timestep outPut} {

    package require pbctools

    # Input:
    set psf $psfPrefix.psf

    switch -- $dim {
    
        "x" {
            set dimension 0
            set a "y"
            set b "z"
        }
        "y" {
            set dimension 1
            set a "x"
            set b "z"
        }
        "z" {
            set dimension 2
            set a "x"
            set b "y"
        }
    }

    # load trajectory and wait until loading completes 
    mol new $psf 
    mol addfile $dcd waitfor all
    set n [ molinfo top get numframes ]
    puts "  Loaded $n frames."
    
    # define first and last frames and stride
    set firstframe 0
    set lastframe $n
    set stride 1
   
    # set output file
    set out [open $outPut a]

    # proceed frame by frame
    for {set i $firstframe} {$i <= $lastframe} {set i [ expr $i + $stride ] } {
        
        set lipid [atomselect top "resname $lipidResname" frame $i]

        set mid [lindex [measure center $lipid] $dimension] 

        set up   [expr $mid + $interval/2.0]
        set down [expr $mid - $interval/2.0]

        set sel_t [atomselect top "name OH2 and ($dim > $down and $dim < $up)" frame $i]
        set sel_1 [atomselect top "name OH2 and ($dim > $down and $dim < $up) and ($a*$a + $b*$b < $R*$R)" frame $i]
        set sel_2 [atomselect top "name OH2 and ($dim > $down and $dim < $up) and ($a*$a + $b*$b > $R*$R)" frame $i]

        puts $out "[$sel_t num]\t[$sel_1 num]\t[$sel_2 num]"
        
        puts "frame $i" 
        #puts "frame $i $mid $up $down $a $b $R [$sel num]" 
        $lipid delete
        $sel_t delete
        $sel_1 delete
        $sel_2 delete
        
    }

    mol delete top

    close $out

}


if {$argc < 8} {
    puts "vmd -dispdev text -e $argv0 -args psfPrefix dcd lipidResname dim interval R timestep outPut"
    exit
}

set psfPrefix    [lindex $argv 0]
set dcd          [lindex $argv 1]
set lipidResname [lindex $argv 2]
set dim          [lindex $argv 3]
set interval     [lindex $argv 4]
set R            [lindex $argv 5]
set timestep     [lindex $argv 6]
set outPut       [lindex $argv 7]

RMSD $psfPrefix $dcd $lipidResname $dim $interval $R $timestep $outPut

exit



