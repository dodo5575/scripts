# jcomer2@uiuc.edu
# Analysis:
set analysisScript pmeAnal.tcl
set outDir potential
set simulFile simulsh.sh

set stride 20

# Input data:
set sim {}
lappend sim {layer_zer lipid ../extra_layer_zer ../output/layer _dopc_zer.dcd {3 4 5} 9984}

lappend sim {layer_q-1 lipid ../anneal_dopc_q-1 ../output/ .dcd {layer2_dopc_q-1} 9984}
lappend sim {layer_q1 lipid ../anneal_dopc_q1 ../output/ .dcd {layer2_dopc_q1} 9984}
lappend sim {layer_q2c lipid ../anneal_dopc_q2c ../output/ .dcd {layer2_dopc_q2c} 9984}

lappend sim {layer_q2s lipid ../anneal_dopc_q2s ../output/layer _dopc_q2s.dcd {1 2} 9984}
lappend sim {layer_q2r lipid ../anneal_dopc_q2r ../output/layer _dopc_q2r.dcd {1 2} 9984}


proc checkExistence {path} {
    if {![file exists $path]} {
	puts "ERROR! The file $path does not exist."
	return 0
    }
    return 1
}

# Check that the files exist.
puts "\nChecking that the necessary files exist..."
set okay 1
set simDcdList {}
foreach s $sim {
    foreach {name moiety structPrefix dcdPrefix dcdSuffix dcdSet dcdFreq} $s {break}    

    checkExistence $structPrefix.psf
    checkExistence $structPrefix.pdb
    #checkExistence $structPrefix.xsc
    set dcdList {}
    foreach d $dcdSet {
	set dcd ${dcdPrefix}${d}${dcdSuffix}
	lappend dcdList $dcd
	if {![checkExistence $dcd]} {set okay 0}
    }
    lappend simDcdList $dcdList
}

if {!$okay} {exit}

set out [open $simulFile w]

# Analyze.
foreach s $sim dcdList $simDcdList {
    foreach {name moiety structPrefix dcdPrefix dcdSuffix dcdSet dcdFreq} $s {break}
    puts "Simulation: $s"
    puts -nonewline $out "vmd -dispdev text -e $analysisScript -args $name $moiety $structPrefix $outDir $dcdFreq $stride"
    foreach dcd $dcdList {
	puts -nonewline $out " $dcd"
    }
    puts $out ""
}

close $out
exit
