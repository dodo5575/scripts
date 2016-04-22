# to use: vmd -dispdev text -e solvate.tcl

package require solvate

solvate deqiang_pore.psf deqiang_pore.pdb -minmax {{-50 -50 -139} {50 50 139}} -o deqiang_sol

file delete combine.psf combine.pdb

quit



