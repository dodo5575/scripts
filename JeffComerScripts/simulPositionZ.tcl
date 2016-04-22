# jcomer2@uiuc.edu
# Analysis:

#set analysisScript analPositionZ.tcl
#set outDir thickness

#set analysisScript analDistZ.tcl
#set outDir thickness

set analysisScript analDistZGap.tcl
set outDir lipid_gap

#set analysisScript analAngle.tcl
#set outDir angle

set simulFile simulsh.sh
set stride 10

# Input data:
set sim {}
lappend sim {layer_pos lipid ../dcd/nw_extra_layer_pos ../dcd/nw_layer _dopc_pos.dcd {0 1 2 3a 4 5 6 7 8 9 10 11 12 13 14 15 16} 9984}
lappend sim {layer_neg lipid ../dcd/nw_extra_layer_neg ../dcd/nw_layer _dopc_neg.dcd {0 1 2 3a 4 5 6 7 8 9 10 11} 9984}
lappend sim {layer_nng lipid ../dcd/nw_extra_layer_nng ../dcd/nw_layer _dopc_nng.dcd {0 1 2 3a 4 5 6 7 8 9 10} 9984}
lappend sim {layer_ngn lipid ../dcd/nw_extra_layer_ngn ../dcd/nw_layer _dopc_ngn.dcd {0 1 2 3 4 5} 9984}
lappend sim {layer_zer lipid ../dcd/nw_extra_layer_zer ../dcd/nw_layer _dopc_zer.dcd {0 1 2 3 4 5} 9984}

#lappend sim {layer_p2e lipid ../dcd/nw_extra_layer_p2e ../dcd/nw_layer _dopc_p2e.dcd {0} 9984}
#lappend sim {layer_p2f lipid ../dcd/nw_extra_layer_p2f ../dcd/nw_layer _dopc_p2f.dcd {0} 9984}

lappend sim {layer_q1 lipid ../dcd/nw_anneal_dopc_q1 ../dcd/nw_layer _dopc_q1.dcd {0 1 2} 9984}
lappend sim {layer_q-1 lipid ../dcd/nw_anneal_dopc_q-1 ../dcd/nw_layer _dopc_q-1.dcd {0 1 2} 9984}
lappend sim {layer_q2c lipid ../dcd/nw_anneal_dopc_q2c ../dcd/nw_layer _dopc_q2c.dcd {0 1 2} 9984}
lappend sim {layer_q0 lipid ../dcd/nw_anneal_dopc_q0 ../dcd/nw_layer _dopc_q0.dcd {0 1} 9984}

lappend sim {layer_aq1 lipid ../dcd/nw_anneal_dopc_q1 ../dcd/nw_layer _dopc_aq1.dcd {0 1} 9984}
lappend sim {layer_aq-1 lipid ../dcd/nw_anneal_dopc_q-1 ../dcd/nw_layer _dopc_aq-1.dcd {0 1} 9984}
lappend sim {layer_aq2c lipid ../dcd/nw_anneal_dopc_q2s ../dcd/nw_layer _dopc_aq2c.dcd {0 1} 9984}
lappend sim {layer_aq0 lipid ../dcd/nw_anneal_dopc_q0 ../dcd/nw_layer _dopc_aq0.dcd {0} 9984}

lappend sim {layer_q2s lipid ../dcd/nw_anneal_dopc_q2s ../dcd/nw_layer _dopc_q2s.dcd {0 1 2} 9984}
lappend sim {layer_q2r lipid ../dcd/nw_anneal_dopc_q2r ../dcd/nw_layer _dopc_q2r.dcd {0 1 2} 9984}


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
