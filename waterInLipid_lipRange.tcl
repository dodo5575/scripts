

proc RMSD {psfPrefix dcd xstPrefix lipidResname dim interval R timestep outPut} {

    package require pbctools

    variable lipR
    variable cutOffFromDNA
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
    mol load psf $psf 
    mol addfile $dcd first 0 last -1 waitfor all
    set n [ molinfo top get numframes ]
    puts "  Loaded $n frames."
    
    # define first and last frames and stride
    set firstframe 0
    set lastframe $n
    set stride 1
   
    # set output file
    set out [open $outPut a]

    # Get the time in nanoseconds.
    set t ""
    set inStream [open ${xstPrefix}.xst r]
    gets $inStream;  gets $inStream;
    while {[gets $inStream line] > 0} {
        lappend t [expr [lindex $line 0] * $timestep * 1.0e-6]
    }
 
    # proceed frame by frame
    for {set i $firstframe} {$i <= $lastframe} {set i [ expr $i + $stride ] } {
        
        set lipid [atomselect top "resname $lipidResname and same residue as (x*x + y*y < $lipR*$lipR)" frame $i]

        set mid [lindex [measure center $lipid] $dimension] 

        set up   [expr $mid + $interval/2.0]
        set down [expr $mid - $interval/2.0]

        set sel_t [atomselect top "name OH2 and ($dim > $down and $dim < $up) and (same residue as within $cutOffFromDNA of nucleic) " frame $i]
        set sel_1 [atomselect top "name OH2 and ($dim > $down and $dim < $up) and (same residue as within $cutOffFromDNA of nucleic)  and ($a*$a + $b*$b < $R*$R)" frame $i]
        set sel_2 [atomselect top "name OH2 and ($dim > $down and $dim < $up) and (same residue as within $cutOffFromDNA of nucleic)  and ($a*$a + $b*$b > $R*$R)" frame $i]

        puts $out "[lindex $t $i]\t[$sel_t num]\t[$sel_1 num]\t[$sel_2 num]"
        
        puts "frame $i" 
        #puts "frame $i $mid $up $down $a $b $R [$sel_t num]" 
        $lipid delete
        $sel_t delete
        $sel_1 delete
        $sel_2 delete
        
    }

    mol delete top

    close $out

}


if {$argc < 8} {
    puts "vmd -dispdev text -e $argv0 -args psfPrefix dcd xstPrefix lipidResname dim interval R timestep outPut"
    exit
}

set psfPrefix    [lindex $argv 0]
set dcd          [lindex $argv 1]
set xstPrefix    [lindex $argv 2]
set lipidResname [lindex $argv 3]
set dim          [lindex $argv 4]
set interval     [lindex $argv 5]
set R            [lindex $argv 6]
set timestep     [lindex $argv 7]
set lipR         [lindex $argv 8]
set cutOffFromDNA         [lindex $argv 9]
set outPut       [lindex $argv 10]

RMSD $psfPrefix $dcd $xstPrefix $lipidResname $dim $interval $R $timestep $outPut

exit



