# Author: Jeff Comer <jcomer2@illinois.edu>

set selText "index 0"
set psf one.psf
set pdb one.pdb
set displayPeriod 100000
# Input:
#set inFile hard_pmf_no.dx
#set inFile diffusion70_pot.dx
set inFile grid0.6.dx 
# Output:
set outFile grid_pos_high.dcd

source $env(HOME)/scripts/vector.tcl
source $env(HOME)/scripts/gridForce.tcl

readDx grid $inFile
puts "Read $inFile, containing [llength $grid(data)] potential values."
puts "Size: $grid(nx) $grid(ny) $grid(nz)"

set molId [mol load psf $psf pdb $pdb]
set sel [atomselect top $selText]

for {set i 0} {$i < $grid(size)} {incr i} {
    animate dup $molId
    set r [indexToWorld grid $i]
    $sel set {x y z} [list $r]
    if {$i % $displayPeriod == 0} { puts "Point $i" }
}

# Delete the first frame (from the pdb).
animate delete beg 0 end 0 $molId

# Write the resulting dcd.
animate write dcd $outFile beg 0 end -1 waitfor all sel $sel $molId

mol delete top
$sel delete
exit
