

proc RMSD {psfPrefix dcdPrefix xstPrefix lipidResname dim R timestep outPut} {

    package require pbctools

    # Input:
    set psf $psfPrefix.psf
    set dcd $dcdPrefix.dcd

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
    mol load psf $psf dcd $dcd
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
    set x ""
    set y ""
    set z ""
    set inStream [open ${xstPrefix}.xst r]
    gets $inStream;  gets $inStream;
    while {[gets $inStream line] > 0} {
        lappend t [expr [lindex $line 0] * $timestep * 1.0e-6]
        lappend x [lindex $line 1]
        lappend y [lindex $line 5]
        lappend z [lindex $line 9]
    }
 
    # proceed frame by frame
    for {set i $firstframe} {$i <= $lastframe} {set i [ expr $i + $stride ] } {
        
        set lipid [atomselect top "resname $lipidResname and name P and ($a*$a + $b*$b > $R*$R)" frame $i]

        set nlipid [$lipid num]

        set area [expr ([lindex $x $i] * [lindex $y $i] - ((atan(1)*4) * $R*$R)) / ($nlipid / 2)]

        puts $out "[lindex $t $i]\t$area"
        
        puts "frame $i" 
        #puts "frame $i $mid $up $down $a $b $R [$sel num]" 
        $lipid delete
        
    }

    mol delete top

    close $out

}


if {$argc < 7} {
    puts "vmd -dispdev text -e $argv0 -args psfPrefix dcdPrefix xstPrefix lipidResname dim R timestep outPut"
    exit
}

set psfPrefix    [lindex $argv 0]
set dcdPrefix    [lindex $argv 1]
set xstPrefix    [lindex $argv 2]
set lipidResname [lindex $argv 3]
set dim          [lindex $argv 4]
set R            [lindex $argv 5]
set timestep     [lindex $argv 6]
set outPut       [lindex $argv 7]

RMSD $psfPrefix $dcdPrefix $xstPrefix $lipidResname $dim $R $timestep $outPut

exit



