# Remove the water from dcd files.
# Author: Jeff Comer <jcomer2@illinois.edu>

set replaceExisting 1
set stride 1
set selText "(not water) and (not resname SIN)"
# Input:
set psf0 pore_none_basepair.psf
set pdb0 pore_none_basepair.pdb
set dcdDir output/
set fileGlob "b*.dcd"
# Output:
set outDir output
set outPrefix nw_
set psf ${outPrefix}${psf0}
set pdb ${outPrefix}${pdb0}

proc trimExtension {name} {
    set ind [string last "." $name]
    return [string range $name 0 [expr {$ind-1}]]
}

proc trimPath {name} {
    set ind [string last "/" $name]
    return [string range $name [expr $ind+1] end]
}

set dcdList [lsort [glob -directory $dcdDir $fileGlob]]
puts "\nFound [llength $dcdList] files.\n"
if {[llength $dcdList] == 0} {
    puts "No match."
    exit
}

if {!$replaceExisting} {
    # Make only the files that don't exist yet.
    set newList {}
    foreach dcd $dcdList {
	if {![file exists $dcd]} { lappend newList $dcd }
    }
    set dcdList $newList
}

# Load the structure and make the selection.
mol load psf $psf0 pdb $pdb0
set sel [atomselect top $selText]
set nAtoms [$sel num]
puts "\nNOTE: Retaining $nAtoms atoms defined by $selText."

set count 0
foreach dcd $dcdList {
   set fileName [trimPath $dcd]
    
    # Skip this file if the output file already exists.
    if {[file exists "$outDir/${outPrefix}${fileName}"]} {
	puts "Skipped $fileName."
	continue
    }

    # Load the trajectory.
    animate delete all
    mol addfile $dcd type dcd step $stride waitfor all
    set nFrames [molinfo top get numframes]
    puts [format "Reading %i frames." $nFrames]
    set last [expr $nFrames - 1]

    # Write the psf once.
    if {$count == 0} {
	$sel writepsf $psf
	$sel writepdb $pdb
    }

    animate write dcd "$outDir/${outPrefix}${fileName}" beg 0 end $last waitfor all sel $sel top
    incr count
}

$sel delete
mol delete top
exit
