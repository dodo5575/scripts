# Author: Jeff Comer <jcomer2@illinois.edu>
source $env(HOME)/scripts/vmdargs.tcl
source $env(HOME)/scripts/useful.tcl

vmdargslist buildStructure.tcl coordFile topFile outDir
set name [trimPath [trimExtension $coordFile]]
set outName $outDir/built_${name}

package require psfgen
topology $topFile
resetpsf

# Load the molecule.
mol new $coordFile
set all [atomselect top all]
set segList [lsort -unique [$all get segname]]
$all delete

# Write the DNA segments.
set tempFileList {}
set segResList {}
foreach seg $segList {
    set sel [atomselect top "segname $seg"]
    set tempFile ${seg}_${name}.pdb
    lappend tempFileList $tempFile

    $sel writepdb $tempFile
    lappend segResList [lsort -unique -integer -index 0 [$sel get {resid resname}]]
    $sel delete
}

# Use psfgen to build the structures.
foreach seg $segList resList $segResList {
    segment $seg {
	#first 5TER
	#last 3TER
	first NONE
	last NONE
	pdb ${seg}_${name}.pdb
    }
    
    foreach res $resList {
	foreach {resid resname} $res {break}
	if {[string equal $resname "THY"] || [string equal $resname "CYT"]} {
	    patch DEO1 $seg:$resid  
	}
	if {[string equal $resname "ADE"] || [string equal $resname "GUA"]} {
	    patch DEO2 $seg:$resid  
	}
    }
    coordpdb ${seg}_${name}.pdb

    # Link to make effectively infinite DNA.
    set resFirst [lindex $resList 0 0]
    set resLast [lindex $resList end 0]
    puts "Linking $seg:$resFirst $seg:$resLast"
    patch LKNA $seg:$resLast $seg:$resFirst 
    regenerate angles dihedrals
}

# Write the completed DNA structures.
guesscoord
writepsf $outName.psf
writepdb $outName.pdb

# Delete the extra files.
foreach f $tempFileList {
    file delete $f
}


exit

