# Author: Jeff Comer <jcomer2@illinois.edu>
set concList {0.1 0.13 1.4 1.9 2 4 6 8}
set volume [expr {40.0*40.0*72.0}]
set prefix bulk_
set eField 0.098281

set stepIons [expr {(100000000*5)/8}]

foreach conc $concList {
    # Convert from molarity to number.
    set num [expr {int(floor($conc*$volume/1.605402e3))}]
    set name c${conc}
    set steps [expr {int(ceil($stepIons/$num))}]

    set fileName ${prefix}${name}.brown
    set out [open $fileName w]

    puts $out "outputName brown"
    puts $out "timestep 2e-5"
    puts $out "steps $steps"
    puts $out "numberFluct 0"
    puts $out "interparticleForce 1"
    puts $out "fullElect 1"
    puts $out "kT 1.0" 
    puts $out "coulombConst 6.156964"
    puts $out "electricField $eField"
    puts $out "outputPeriod 2000"
    puts $out "outputFormat  dcd"
    puts $out "cutoff 40.0"
   
    puts $out ""
    puts $out "particle POT"
    puts $out "num $num"
    puts $out "gridFile zero_pmf.dx"
    puts $out "diffusion 227"
    puts $out "charge 1"
    puts $out "radius 1.76375"
    puts $out "eps 0.0870"

    puts $out ""
    puts $out "particle CLA"
    puts $out "num $num"
    puts $out "gridFile zero_pmf.dx"
    puts $out "diffusion 241"
    puts $out "charge -1"
    puts $out "radius 2.27"
    puts $out "eps 0.150"

    puts $out ""
    puts $out "tabulatedPotential 1"
    puts $out "tabulatedFile 0@0@long_pot-pot.dat"
    puts $out "tabulatedFile 0@1@long_pot-chl.dat"
    puts $out "tabulatedFile 1@1@long_chl-chl.dat"
}
