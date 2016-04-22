# averageStructure.tcl --- produce an averaged structure

source $env(SCRIPTS)/Procs.tcl

### PARAMETERS ###
#
set aligntext "protein name CA"


### ARGUMENTS ###
#
vmdargs averageStucture.tcl psf dcd first last outbase


### MAIN ###
#
mol load psf $psf
mol addfile $dcd first $first last $last waitfor all
set nframes [molinfo top get numframes]

set coords_sum {}
set Lz_sum 0.0

set alignsel [atomselect top $aligntext]
set refsel [atomselect top $aligntext frame 0]
set all [atomselect top all]

for { set frame 0 } { $frame < $nframes } { incr frame } {
    progressbar $frame $nframes
    
    animate goto $frame
    $refsel frame 0
#     $alignsel frame $frame
#     $all frame $frame
    
    set mat [measure fit $alignsel $refsel]
    $all move $mat
    
    set coords_sum [vecadd_list $coords_sum [$all get {x y z}]]
    set Lz_sum [expr $Lz_sum + [molinfo top get c]]
}
set coords [vecscale_list $coords_sum [expr {1.0/$nframes}]]; silent
set Lz [expr {$Lz_sum/$nframes}]; silent

puts "<Lz> = $Lz"
set ouch [open $outbase.pdb.Lz w]
puts $ouch "<Lz> = $Lz"
close $ouch

$all set {x y z} $coords

$all writepdb $outbase.pdb

exit
