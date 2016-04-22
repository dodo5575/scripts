# Get the system size using measure minmax.
# to use: vmd -dispdev text -e getSize.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set x0 20
set x1 35
set y0 -15
set y1 35
set z0 -30
set z1 30
set selText "resname SIO2"

# Delete all representations.
set n [molinfo top get numreps]
for {set i 0} {$i < $n} {incr i} { mol delrep 0 top }

set dalton 1.66054e-27; # kg

set selStuff "($selText) and x > $x0 and x < $x1 and y > $y0 and y < $y1 and z > $z0 and z < $z1"
set sel [atomselect top $selStuff]
puts "Number: [$sel num]"
set volume [expr ($x1-$x0)*($y1-$y0)*($z1-$z0)]
puts "Volume: $volume A^3"
set mass [measure sumweights $sel weight mass]
puts "Mass: $mass Da"
puts "Density: [expr $mass/$volume] Da/A^3"
puts "Density: [expr 1e3*$mass*$dalton/($volume*1e-24)] g/cm^3"

# Make a representation.
mol representation VDW 0.800000 15.000000
mol color ColorID 0
mol selection $selStuff
mol material Opaque
mol addrep top

mol representation Points 1.0
mol color ColorID 1
mol selection $selText
mol material Opaque
mol addrep top

