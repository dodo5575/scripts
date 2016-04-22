#!/usr/bin/tclsh
# Split a file of mini-trajectories into many files differing by position in space. 
# Author: Jeff Comer <jcomer2@illinois.edu>

source vector.tcl
source gridForce.tcl

if {$argc != 3} {
    puts "Usage: ./assignDiffusionToGrid.tcl indexFile gridFile outName"
    exit
}

set poreZ0 12
set poreS0 10
set farZ0 20

set indexFile [lindex $argv 0]
set gridFile [lindex $argv 1]
set outName [lindex $argv 2]
set kT 0.5862292
# Convert diffusivities to A^2/ns. Was in A^2/fs.
set convFactor 1e6

proc trimExtension {name} {
    set ind [string last "." $name]
    return [string range $name 0 [expr {$ind-1}]]
}

proc trimPath {name} {
    set ind [string last "/" $name]
    return [string range $name [expr $ind+1] end]
}

proc diffusionFromVel {velList dt} {
    # Integrate the correlation function to obtain the diffusivity.
    set sum 0.0
    foreach corr $velList {
	set sum [expr $sum + $corr*$dt]
    }
    
    # Scale by the number of dimensions.
    set diffuse [expr {$sum/3.0}]

    return $diffuse
}

proc diffusionFromForce {forceList dt kT} {
    # Integrate the correlation function to obtain the diffusivity.
    set sum 0.0
    foreach corr $forceList {
	set sum [expr $sum + $corr*$dt]
    }
   
    # Compute the diffusivity.
    if {$sum == 0.0} {
	return -1
    }
    set diffuse [expr {$kT*$kT/$sum}]

    return $diffuse
}

proc writeCorrFunction {fileName l dt} {
    set out [open $fileName w]
    set step 0
    foreach i $l {
	puts $out "[expr {$dt*$step}] $i"
	incr step
    }
    close $out
    return
}

proc readFunction {inFile} {
    # Open the files.
    set in [open $inFile r]
    set timList {}
    set posList {}
    while {[gets $in line] >= 0} {
	if {[string match "#*" $line]} { continue }

	set item [concat $line]
	if {[llength $item] != 2} { 
	    puts "Warning: Invalid line `$line'."
	    continue 
	}
	lappend timList [lindex $item 0]
	lappend posList [lindex $item 1]
    }
    close $in

    return [list $timList $posList]
}

proc meanTrajList {listList} {
    set meanList {}
    set n [llength [lindex $listList 0]]
    
    for {set i 0} {$i < $n} {incr i} {
	set sum 0.0
	set count [llength $listList]
	foreach item $listList {
	    set sum [expr {$sum + [lindex $item $i]}]
	}
	lappend meanList [expr $sum/$count]
    }
    
    return $meanList
}

proc meanWeightTrajList {listList weightList} {
    set meanList {}
    set n [llength [lindex $listList 0]]
    
    for {set i 0} {$i < $n} {incr i} {
	set sum 0.0
	set weightSum 0.0
	foreach item $listList w $weightList {
	    set sum [expr {$sum + $w*[lindex $item $i]}]
	    set weightSum [expr {$weightSum + $w}]
	}
	lappend meanList [expr $sum/$weightSum]
    }
    
    return $meanList
}

# Get the grid geometry.
readDx grid $gridFile
copyGridDim grid diffFGrid
copyGridDim grid diffVGrid
puts "Read the grid with dimensions $grid(nx) $grid(ny) $grid(nz)."

# Lists for the special diffusivity regions.
set poreWeight {}
set poreCorrV {}
set poreCorrF {}
set farWeight {}
set farCorrV {}
set farCorrF {}

# Read the index file.
set in [open $indexFile]
while {[gets $in line] >= 0} {
    if {[string length $line] <= 1} { continue }
    
    if {[string match "#*" $line]} { continue }

    set item [concat $line]
    if {[llength $item] < 2} {
	puts "Warning: Invalid line `$line'" 
	continue 
    }

    # Get the info from the index file.
    set filePrefix [trimExtension [lindex $item 0]]
    set weight [lindex $item 1]

    # Read the correlation functions.    
    set stuff [readFunction $filePrefix.vcorr]
    set timList [lindex $stuff 0]
    set velList [lindex $stuff 1]
    set stuff [readFunction $filePrefix.fcorr]
    set forceList [lindex $stuff 1]
    set dt [expr {double([lindex $timList 1]-[lindex $timList 0])}]

    # Find the position of the node in the grid.
    set index [regexp -inline {[0-9]*$} $filePrefix]
    set pos [indexToWorld grid $index]
    puts "pos $filePrefix $index $pos"
    foreach {x y z} $pos { break }
    set s [expr {sqrt($x*$x + $y*$y)}]
    
    # Compute the diffusivity of this grid node.
    set diffV [expr {$convFactor*[diffusionFromVel $velList $dt]}]
    set diffF [expr {$convFactor*[diffusionFromForce $velList $dt $kT]}]
    puts "index $index diffV $diffV diffF $diffF pos $x $y $z"

    # Add the diffusivity to the grids.
    lset diffVGrid(data) $index $diffV
    lset diffFGrid(data) $index $diffF

    # Accumulate the autocorrelation functions for the special regions.
    if {$s < $poreS0 && abs($z) < $poreZ0} {
	lappend poreWeight $weight
	lappend poreCorrV $velList
	lappend poreCorrF $forceList
    }
    if {abs($z) > $farZ0} {
	lappend farWeight $weight
	lappend farCorrV $velList
	lappend farCorrF $forceList
    }
    
}
close $in

# Write the resulting grids.
set outFileV ${outName}_diffVel.dx
writeDx diffVGrid ${outName}_diffVel.dx
puts "Wrote the diffusivity grid derived from the velocity to `$outFileV'." 

set outFileF ${outName}_diffForce.dx
writeDx diffFGrid ${outName}_diffForce.dx
puts "Wrote the diffusivity grid derived from the force to `$outFileF'." 

# Compute the special region diffusivities.
set poreMeanCorrV [meanWeightTrajList $poreCorrV $poreWeight]
set poreDiffV [expr {$convFactor*[diffusionFromVel $poreMeanCorrV $dt]}]
puts "poreDiffV $poreDiffV"
writeCorrFunction ${outName}_poreCorrV.dat $poreMeanCorrV $dt

set poreMeanCorrF [meanWeightTrajList $poreCorrF $poreWeight]
set poreDiffF [expr {$convFactor*[diffusionFromForce $poreMeanCorrF $dt $kT]}]
puts "poreDiffF $poreDiffF"
writeCorrFunction ${outName}_poreCorrF.dat $poreMeanCorrF $dt

set farMeanCorrV [meanWeightTrajList $farCorrV $farWeight]
set farDiffV [expr {$convFactor*[diffusionFromVel $farMeanCorrV $dt]}]
puts "farDiffV $farDiffV"
writeCorrFunction ${outName}_farCorrV.dat $farMeanCorrV $dt

set farMeanCorrF [meanWeightTrajList $farCorrF $farWeight]
set farDiffF [expr {$convFactor*[diffusionFromForce $farMeanCorrF $dt $kT]}]
puts "farDiffF $farDiffF"
writeCorrFunction ${outName}_farCorrF.dat $farMeanCorrF $dt

exit
