#!/usr/bin/tclsh
# jcomer2@illinois.edu
if {[llength $argv] != 3} {
    puts "$argv0 diffNodeFile structPrefix outPrefix"
    exit
}
set diffNodeFile [lindex $argv 0]
set structPrefix [lindex $argv 1]
set outPrefix [lindex $argv end]

set selText "nucleic and noh"
set segName M
set name C
set resName M

# Construct a pdb line from everything.
proc makePdbLineBeta {index segName resId name resName r beta} {
    set template "ATOM     42  HN  ASN X   4     -41.083  17.391  50.684  0.00  0.00      P1    "

    foreach {x y z} $r {break}
    set record "ATOM  "
    set si [string range [format "     %5i " $index] end-5 end]
    if {[string length $name] < 4} {
	set name [string range " $name    " 0 3]
    } else {
	set name [string range $name 0 3]
    }
    set resName [string range " $resName    " 0 3]
    set temp0 " [string index $segName 0]"
    set resId [string range "    $resId"  end-3 end]
    set temp1 [string range $template  26 29]
    set sx [string range [format "       %8.3f" $x] end-7 end]
    set sy [string range [format "       %8.3f" $y] end-7 end]
    set sz [string range [format "       %8.3f" $z] end-7 end]
    set temp2 [string range $template 54 59]
    set beta [string range [format "       %6.2f" $beta] end-5 end]
    set temp3 [string range $template 66 71]
    set segName [string range "$segName    "  0 3]
    set tempEnd [string range $template 76 end]

    # Construct the pdb line.
    return "${record}${si}${name}${resName}${temp0}${resId}${temp1}${sx}${sy}${sz}${temp2}${beta}${temp3}${segName}${tempEnd}"
}

# Just read a space delimited data file.
proc readData {fileName} {
    set in [open $fileName r]
    
    set r {}
    while {[gets $in line] >= 0} {
	if {[string match "#*" $line]} {continue}
	if {[string length $line] < 2} {continue}

	set tok [concat $line]
	lappend r $tok
    }

    close $in
    return $r
}

proc trimPath {name} {
    set ind [string last "/" $name]
    return [string range $name [expr $ind+1] end]
}

proc extractPath {name} {
    set ind [string last "/" $name]
    return [string range $name 0 $ind]
}

proc quiet {} {
}

# Get the window data and other stuff.
set diffuseData [readData $diffNodeFile]; quiet

# Get the distance between window nodes and DNA.
mol load psf $structPrefix.psf pdb $structPrefix.pdb
set sel [atomselect top $selText]
set atomPosList [$sel get {x y z}]; quiet
set n [$sel num]
$sel delete
mol delete top

set distList {}
foreach item $diffuseData {
    set meanPos [lrange $item 1 3]

    # Find the nearest atom to this window.
    set pos [lindex $atomPosList 0]
    set d [veclength [vecsub $pos $meanPos]]

    puts $d
    set minDist $d
    foreach pos $atomPosList {
	set d [veclength [vecsub $pos $meanPos]]
	if {$d < $minDist} { 
	    set minDist $d
	}
    }
    
    lappend distList $minDist
}

# Write the output.
set outDir [extractPath $outPrefix]
set outName [trimPath $outPrefix]

# Sort the windows.
set elSort {}
set out [open $outDir/$outName.pdb w]
set i 1
foreach item $diffuseData dist $distList {
    set r [lrange $item 1 3]
    set d [lrange $item 4 6]
    set dMean [expr {([lindex $d 0] + [lindex $d 1] + [lindex $d 2])/3.0}]

    # Get the minimum and maximum component.
    set dMin [lindex $d 0]
    set dMax [lindex $d 0]
    foreach dComp $d {
	if {$dComp < $dMin} { set dMin $dComp }
	if {$dComp > $dMax} { set dMax $dComp }
    }
    set errMin [expr {$dMean-$dMin}]
    set errMax [expr {$dMax-$dMean}]

    set dz [lindex $d 2]
    set name "C[expr {int(floor($dMean))}]"

    lappend meanSort [list $dist $dMean $errMin $errMax]
    lappend elSort [concat $dist $d]

    set lin [makePdbLineBeta $i $segName $i $name $resName $r $dMean]
    puts $out $lin
    incr i
}
close $out

set meanSort [lsort -real -index 0 $meanSort]; quiet
set out [open $outDir/profileMean_${outName}.dat w]; quiet
foreach item $meanSort {
    puts $out "$item"
}
close $out

set elSort [lsort -real -index 0 $elSort]; quiet

set outX [open $outDir/profileDx_${outName}.dat w]; quiet
set outY [open $outDir/profileDy_${outName}.dat w]; quiet
set outZ [open $outDir/profileDz_${outName}.dat w]; quiet
foreach item $elSort {
    puts $outX "[lindex $item 0] [lindex $item 1]"
    puts $outY "[lindex $item 0] [lindex $item 2]"
    puts $outZ "[lindex $item 0] [lindex $item 3]"
}
close $outX
close $outY
close $outZ

exit
