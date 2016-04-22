# to use: vmd -dispdev text -e solvate.tcl

package require solvate

set sys s-dna
#solvate pore6_${sys}.psf pore6_${sys}.pdb -minmax {{-80 -80 -150} {80 80 150}} -s W -spsf water/water.psf -spdb water/water.pdb -stop water/tip3p.top -ks "name OH2" -ws 40 -b 1.7 -o pore6_${sys}_sol
solvate pore6_${sys}.psf pore6_${sys}.pdb -minmax {{-80 -80 -170} {80 80 170}} -s W -ks "name OH2" -b 1.7 -o pore6_${sys}_sol

file delete combine.psf combine.pdb

quit

