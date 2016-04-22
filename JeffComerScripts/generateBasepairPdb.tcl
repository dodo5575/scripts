set inPdb approach_sep0.txt.pdb
set outPdb reorder.pdb
set nBasepairs 5

source $env(HOME)/scripts/useful.tcl
mol load pdb $inPdb
set out [open $outPdb w]

set ind 1
for {set i 0} {$i < $nBasepairs} {incr i} {
    set seg D${i}
    set segName ADNA

    set sel [atomselect top "segname $seg and resid 1"]
    set posList [$sel get {x y z}]
    set nameList [$sel get name]
    set resNameList [$sel get resname]

    foreach pos $posList name $nameList resName $resNameList {
	puts $out [makePdbLineBeta $ind $segName [expr {$i+1}] $name $resName $pos 1.0]
	incr ind
    }
    
}

for {set i 0} {$i < $nBasepairs} {incr i} {
    set seg D[expr {$nBasepairs-$i-1}]
    set segName BDNA

    set sel [atomselect top "segname $seg and resid 2"]
    set posList [$sel get {x y z}]
    set nameList [$sel get name]
    set resNameList [$sel get resname]

    foreach pos $posList name $nameList resName $resNameList {
	puts $out [makePdbLineBeta $ind $segName [expr {$i+1}] $name $resName $pos 1.0]
	incr ind
    }
    
}

close $out
exit
