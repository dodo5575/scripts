# Author: Jeff Comer <jcomer2@illinois.edu>
set nameList {at_neg at_pos gc_neg gc_pos gcm_neg gcm_pos no_neg no_pos}
set prefix best1_
set gridPrefix best
set numberFluct 0
set electricField 0.098281
set sep ""
set steps 600000000

foreach name $nameList {
    set fileName ${prefix}${name}.brown

    set out [open $fileName w]

    if {[string match "*neg*" $name]} {
	set eField [expr -$electricField]
    } else {
	set eField $electricField
    }
    
    if {[string match "*at_*" $name]} {
	set pair at
    } elseif {[string match "gc_*" $name]} {
	set pair gc
    } elseif {[string match "*gcm_*" $name]} {
	set pair gcm
    } else {
	set pair no
    }

    if {[string equal $pair no]} {
	set gridPot  ${gridPrefix}_pmf_${pair}.dx
	set gridChl  ${gridPrefix}_pmf_${pair}.dx
    } else {
	set gridPot  ${gridPrefix}_pmf_${pair}_pot.dx
	set gridChl  ${gridPrefix}_pmf_${pair}_chl.dx
    }

    puts $out "outputName $sep ${prefix}${name}"
    puts $out "timestep $sep 2e-5"
    puts $out "steps $sep $steps"
    puts $out "numberFluct $sep $numberFluct"
    puts $out "interparticleForce $sep 1"
    puts $out "fullElect $sep 1"
    puts $out "kT $sep 1.0" 
    puts $out "coulombConst $sep 6.156964"
    puts $out "electricField $sep $eField"
    puts $out "outputPeriod $sep 2000"
    puts $out "outputPdb $sep 1"
    puts $out "cutoff $sep 10.0"

    puts $out ""
    puts $out "particle $sep POT"
    puts $out "num $sep 7"
    puts $out "gridFile $sep $gridPot"
    puts $out "diffusion $sep 227"
    puts $out "charge $sep 1"
    puts $out "radius $sep 1.76375"
    puts $out "eps $sep 0.0870" 

    puts $out ""
    puts $out "particle $sep CLA"
    puts $out "num $sep 5"
    puts $out "gridFile $sep $gridChl"
    puts $out "diffusion $sep 241"
    puts $out "charge $sep -1"
    puts $out "radius $sep 2.27"
    puts $out "eps $sep 0.150"
}
