# Author: Jeff Comer <jcomer2@illinois.edu>

# Dave's stuff.
# ~david/Work/Scripts/Procs.tcl
proc vmdargs { script args } {
    upvar argv argv_local
    
    if { [llength $argv_local] != [llength $args] } {
	puts -nonewline "Usage: vmd -dispdev text -e $script -args"
	foreach arg $args {
	    puts -nonewline " $arg"
	}
	puts ""
	exit -1
    }
    
    foreach arg $args {
	upvar $arg arg_local
	set arg_local [shift argv_local]
    }
}

proc shift { stack } {
    upvar $stack list
    set value [lindex $list 0]
    set list [lrange $list 1 end]
    return $value
}


vmdargs gridToDcdBig.tcl name

set selText "ions"
set psf $name.psf
set pdb $name.pdb
set displayPeriod 20000
set maxDcdFrames 100000
# Input:
#set inFile hard_pmf_no.dx
#set inFile diffusion70_pot.dx
set inFile grid0.3.dx 
# Output:
set outPrefix dcd/grid_${name}

source $env(HOME)/scripts/vector.tcl
source $env(HOME)/scripts/gridForce.tcl

readDx grid $inFile
puts "Read $inFile, containing [llength $grid(data)] potential values."
puts "Size: $grid(nx) $grid(ny) $grid(nz)"

set molId [mol load psf $psf pdb $pdb]
set sel [atomselect top $selText]
set all [atomselect top all]
set pos0 [$all get {x y z}]

# Delete the first frame (the pdb).
animate delete beg 0 end 0 $molId
	
set count 0
set nDcd 0
for {set i 0} {$i < $grid(size)} {incr i} {
    animate dup $molId
    
    # Set the new position.
    $all set {x y z} $pos0
    set r [indexToWorld grid $i]
    $sel set {x y z} [list $r]

    if {$i % $displayPeriod == 0} { puts "Point $i" }

    incr count
    if {$count >= $maxDcdFrames} {	
	set outFile $outPrefix.$nDcd.dcd
	incr nDcd

	# Write the resulting dcd.
	animate write dcd $outFile beg 0 end -1 waitfor all sel $all $molId
	animate delete all
	set count 0
	puts "Next file!"
    }
}

# Write whatever is left.
puts "Last file!"
set outFile $outPrefix.$nDcd.dcd
animate write dcd $outFile beg 0 end -1 waitfor all sel $all $molId

$all delete

mol delete top
$sel delete
