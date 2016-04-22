# Author: Jeff Comer <jcomer2@illinois.edu>
set topFile muthuDna.top
set n 80
set l0 25; # in A
set segName DNA
set name D
set resName DNA
set dl 2; # in A
# Ouptut:
set outName wiggly80
set tempName ${outName}a

# Construct a pdb line from everything.
proc makePdbLineFull {index segName resId name resName r} {
    foreach {x y z} $r {break}
    set record "ATOM  "
    set si [string range [format "     %5i " $index] end-5 end]
    if {[string length $name] < 4} {
        set name [string range " $name    " 0 3]
    } else {
	set name [string range $name 0 3]
    }
    set name [string range $name 0 3]
    set resName [string range " $resName    " 0 3]
    set temp0 " [string index $segName 0]"
    set resId [string range "    $resId"  end-3 end]
    set temp1 "    "
    set sx [string range [format "       %8.3f" $x] end-7 end]
    set sy [string range [format "       %8.3f" $y] end-7 end]
    set sz [string range [format "       %8.3f" $z] end-7 end]
    set temp2 "  1.00  1.00      "
    set segName [string range "$segName    "  0 3]
    set tempEnd "   "

    # Construct the pdb line.
    return "${record}${si}${name}${resName}${temp0}${resId}${temp1}${sx}${sy}${sz}${temp2}${segName}${tempEnd}"
}

# Open the file.
set out [open $tempName.pdb w]

set x0 0.0
set y0 0.0
set z0 [expr -0.5*$l0*$n]

# Write the beads.
for {set i 1} {$i <= $n} {incr i} {
    set dx [expr 2.0*(rand()-0.5)*$dl]
    set dy [expr 2.0*(rand()-0.5)*$dl]
    set r [list [expr $x0+$dx] [expr $y0+$dy] [expr $z0 + $i*$l0]]
    set lin [makePdbLineFull $i $segName $i $name $resName $r]
    puts $out $lin
}

puts $out END
close $out

# Have psfgen build the structure.
package require psfgen
topology $topFile 
resetpsf
#psfalias

segment $segName {
    auto angles
    first none
    last none
    pdb $tempName.pdb
}
coordpdb $tempName.pdb   
#guesscoord

writepsf $outName.psf
writepdb $outName.pdb
exit
