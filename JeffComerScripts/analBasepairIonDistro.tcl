# Author: jcomer2@illinois.edu
# output is in mol/L

if {$argc < 7} {
    puts "$argv0 name moiety structPrefix outputDir dcdFreq stride dcdFile0 [dcdFile1...]"
    exit
}
set name [lindex $argv 0]
set moiety [lindex $argv 1]
set structPrefix [lindex $argv 2]
set outDir [lindex $argv 3]
set dcdFreq [lindex $argv 4]
set stride [lindex $argv 5]
set dcdList [lrange $argv 6 end]

source $env(HOME)/scripts/vector.tcl
source $env(HOME)/scripts/gridForce.tcl

proc compute {name moiety structPrefix dcdList dcdFreq outDir stride} {
    # Parameters:
    set ion $moiety
    set ionText "name ${ion}"
    set localGrid basepair_local1.dx
    set segA ADNA
    set resA 1
    set segB BDNA
    set resB 111
    set nBasepairs 20

    set startFrame 0
    set displayPeriod 50
    set timestep 1.0
    if {$stride < 1} {set stride 1}

    set avogadro 6.0221367e23
    set convert [expr {1e27/$avogadro}]
    
    # Input:
    set pdb [lindex $dcdList 0].pdb
    set xsc $structPrefix.xsc
    set dnaPsf $structPrefix.psf
    set dnaPdb $structPrefix.pdb    

    # Read the system size from the xsc file.
    set in [open $xsc r]
    foreach line [split [read $in] "\n"] {
	if {![string match "#*" $line]} {
	    set param [split $line]
	    puts $param
	    set ex [lrange $param 1 3]
	    set ey [lrange $param 4 6]
	    set ez [lrange $param 7 9]
	    break
	}
    }
    close $in

    # Load the grids.
    newGridBox sys $ex $ey $ez 5.0
    readDx count $localGrid
    copyGridDim count density
  
    puts "Grid dimensions: $count(nx) $count(ny) $count(nz)"
    puts "$count(delta)"
  
    puts "$timestep $dcdFreq $stride"
    # Get the time change between frames in nanoseconds.
    set dt [expr {1.0e-6*$timestep*$dcdFreq*$stride}]

    # Load the DNA system.
    set dnaMol [mol load psf $dnaPsf pdb $dnaPdb]
    
    # Load the ion system.
    set sysMol [mol load pdb $pdb]
    set ionSel [atomselect $sysMol $ionText]

    # Make a list of the residues.
    set resAList {}
    set segAList {}
    set resBList {}
    set segBList {}
    for {set i 0} {$i < $nBasepairs} {incr i} {
	set a [expr {$resA + $i}]
	set b [expr {$resB - $i}]
	lappend segAList $segA
	lappend resAList $a
	lappend segBList $segB
	lappend resBList $b

	set selA [atomselect $dnaMol "segname $segA and resid $a"]
	set selB [atomselect $dnaMol "segname $segB and resid $b"]
	if {[$selA num] < 3} {
	    puts "ERROR! Residue $segA:$a does not exist."
	    exit
	}
	if {[$selB num] < 3} {
	    puts "ERROR! Residue $segB:$b does not exist."
	    exit
	}
	$selA delete
	$selB delete
    }
        
    # Loop over the dcd files.
    set nFrames0 0
    set frameCount 0
    foreach dcdFile $dcdList {
	# Load the trajectory.
	animate delete all
	mol addfile $dcdFile type dcd step $stride waitfor all
	set nFrames [molinfo $sysMol get numframes]
	puts [format "Reading %i frames." $nFrames]

	# Move forward, computing at each step.
	set pos0 {}
	for {set f $startFrame} {$f < $nFrames} {incr f} {
	    molinfo $sysMol set frame $f

	    set posList [$ionSel get {x y z}]

	    # Run through the basepairs.
	    foreach segA $segAList resA $resAList segB $segBList resB $resBList {
		#Get the transformation to the local basis.
		set basis [getBasepairBasisBest $segA $resA $segB $resB $dnaMol]
		set cen [getBasepairPos $segA $resA $segB $resB $dnaMol]
		set basisInv [matTranspose $basis]

		foreach pos $posList {
		    # Wrap nearest to the center.
		    set p [wrapDiff sys [vecsub $pos $cen]]
		    # Rotate.
		    set posLoc [vecTransform $basisInv $p]

		    # Find the closest local grid point.
		    # Increment the count.
		    if {[inGrid count $posLoc]} {
			set j [nearestIndex count $posLoc]
			set c0 [lindex $count(data) $j]
			lset count(data) $j [expr {$c0 + 1.0}]
		    }
		}
		incr frameCount
	    }
	    
	    # Update the display.
	    if {$f % $displayPeriod == 0} {
		puts [format "FRAME %i" $f ]
	    }
	    set pos0 $pos
	}
	set nFrames0 [expr $nFrames+$nFrames0]
    }

    mol delete $dnaMol
    mol delete $sysMol

    # Compute the mean number density.
    set volume [getVolume count]
    set meanDensity {}
    set meanCount {}
    foreach p $count(data) {
	lappend meanCount [expr {$p/$frameCount}]
	lappend meanDensity [expr {$convert*$p/$frameCount/$volume}]
    }
    set count(data) $meanCount
    set density(data) $meanDensity

    # Write the results.
    writeDx count $outDir/number${ion}_${name}.dx
    writeDx density $outDir/density${ion}_${name}.dx
}


# Get the position of an atom.
proc getPos {segName resId name {mole top}} {
    set sel [atomselect $mole "segname $segName and resid $resId and name $name"]
    set n [$sel num]
    
    if {$n < 1} {
	puts "Warning! Atom ${segName}:${resId}:${name} does not exist."
	return [list 0.0 0.0 0.0]
    } elseif {$n > 1} {
	puts "Warning! Atom ${segName}:${resId}:${name} in not unique."
    }

    set r [lindex [$sel get {x y z}] 0]
    $sel delete
    return $r
}

# Define a basis for a base.
proc getBaseNormal {segName resId {mole top}} {
    set selText "segname $segName and resid $resId"
    # Get the residue name.
    set sel [atomselect $mole $selText]
    set resName [lindex [$sel get resname] 0]
    $sel delete

    # Get the hexagon ring basis.
    if {[string equal ADE $resName] || [string equal GUA $resName]} {
	set selX0 [atomselect $mole "($selText) and name C4"]
	set selX1 [atomselect $mole "($selText) and name N1"]
	set selY0 [atomselect $mole "($selText) and name N3"]
	set selY1 [atomselect $mole "($selText) and name C5"]
    } else {
	set selX0 [atomselect $mole "($selText) and name N1"]
	set selX1 [atomselect $mole "($selText) and name C4"]
	set selY0 [atomselect $mole "($selText) and name C2"]
	set selY1 [atomselect $mole "($selText) and name C6"]
    }

    set rX0 [lindex [$selX0 get {x y z}] 0]
    set rX1 [lindex [$selX1 get {x y z}] 0]
    set rY0 [lindex [$selY0 get {x y z}] 0]
    set rY1 [lindex [$selY1 get {x y z}] 0]

    $selX0 delete
    $selX1 delete
    $selY0 delete
    $selY1 delete

    set ex [vecsub $rX1 $rX0]
    set ex [vecscale [expr 1.0/[veclength $ex]] $ex]
    set ey [vecsub $rY1 $rY0]
    set ey [vecsub $ey [vecscale [vecdot $ey $ex] $ex]]
    set ey [vecscale [expr 1.0/[veclength $ey]] $ey]
    set ez [veccross $ex $ey]

    return $ez
}

# Define a basis for a base.
proc getBasepairDirection {segA resA segB resB {mole top}} {
    set car1A [getPos $segA $resA "C1'" $mole]
    set car1B [getPos $segB $resB "C1'" $mole]

    return [vecUnit [vecsub $car1B $car1A]]
}

proc getBasepairBasisBest {segA resA segB resB {mole top}} {
    set zA [getBaseNormal $segA $resA $mole]
    set zB [getBaseNormal $segB $resB $mole]
    set ez [vecUnit [vecsub $zA $zB]]

    set ex [getBasepairDirection $segA $resA $segB $resB $mole]
    set ex [vecUnit [vecsub $ex [vecscale [vecdot $ex $ez] $ez]]]
    
    set ey [veccross $ez $ex]
    return [matTranspose [list $ex $ey $ez]]
}

# Get the standard position of a DNA basepair.
proc getBasepairPos {segA resA segB resB {mole top}} {
    set car1A [getPos $segA $resA "C1'" $mole]
    set car1B [getPos $segB $resB "C1'" $mole]

    return [vecscale 0.5 [vecadd $car1A $car1B]] 
}

proc splitLast {name char} {
    set ind [string last $char $name]
    return [list [string range $name 0 [expr {$ind-1}]] [string range $name $ind end]]
}

compute $name $moiety $structPrefix $dcdList $dcdFreq $outDir $stride
exit
