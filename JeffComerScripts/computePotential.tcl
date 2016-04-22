# Author: Jeff Comer <jcomer2@illinois.edu>
# Analysis:

set analysisScript electrodeTraj.tcl
set outDir electrode_good
set probeName bottom
set probePdb ${probeName}Probe.pdb

#Input data:
set sim {}
lappend sim {eco_grisha_orig_2V_100ps ecoRI ../ecoRI/grisha/nw_1CKQ_deep_in_pore27_0.1M ../ecoRI/grisha/ {""} 100000}

lappend sim {eco_grisha_2Vb ecoRI ../ecoRI/2V_cont/nw_ecoRI_2V ../ecoRI/2V_cont/dcd {0.0 1.0 1.1 1.2 1.3 2.0 3.0 3.1 3.2 3.3 4.0 4.1 4.2 4.3 5.0 5.1 6.0 6.1 7.0 7.1 7.2 8.0 8.1 8.2} 5000}
lappend sim {eco_start_2.5V ecoRI ../ecoRI/start/nw_eco_start_fold ../ecoRI/start/dcd {0.0 0.1 1.0 1.1 2.0 2.1 3.0 3.1 4.0 4.1 5.0 5.1 6.0 6.1 7.0 7.1} 5000}
lappend sim {eco_start_4V ecoRI ../ecoRI/start/nw_eco_start_fold ../ecoRI/start/dcd {0.0 1.0 1.1 1.2} 5000}

lappend sim {bam_nonspec_0.5V bamHI ../bamHI/nonspecific/nw_bam_eq5_ion2 ../bamHI/nonspecific/dcd {0.0 0.1 1.0 2.0 2.1 2.2 3.0 5.0 5.1 5.2 6.0 6.1 6.2} 5000}
lappend sim {bam_nonspec_1Va bamHI ../bamHI/nonspecific/nw_bam_eq5_ion2 ../bamHI/nonspecific/dcd {0.0 0.1 1.0 1.1 2.0 2.1 3a.0 4.0 4.1 5.0 5.1} 5000}
lappend sim {bam_nonspec_1.5V bamHI ../bamHI/nonspecific/nw_bam_eq5_ion2 ../bamHI/nonspecific/dcd {0.0 0.1 1.0 1.1 2.0 2.1} 5000}
lappend sim {bam_2Va bamHI ../bamHI/nonspecific/nw_bam_all ../bamHI/nonspecific/dcd {0.0 1.0 1.1 1.2 2.0 2.1} 5000}
lappend sim {bam_specific_1.5V bamHI ../bamHI/specific/nw_bam_spec_sys ../bamHI/specific/dcd {0.0 0.1 1.0 1.1 2.0 3.0 3.1 3.2 4.0 4.1 5.0 5.1 5.2 6.0 6.1} 5000}
lappend sim {bam_specific_2V bamHI ../bamHI/specific/nw_bam_spec_sys ../bamHI/specific/dcd {0.0 1.0 1.1 2.0 2.1 3.0 3.1} 5000}
lappend sim {bam_spec_nogrid_2V bamHI ../bamHI/specific/nw_bam_spec_sys ../bamHI/specific/dcd {0.0 0.1 0.2 1.0 2a.0 2a.1 3.0 3.1 3.2 4.0 4.1 4.2} 5000}

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
foreach s $sim {
    foreach {name moiety structPrefix dir dcdSet dcdFreq} $s {break}    

    checkExistence $structPrefix.psf
    checkExistence $structPrefix.pdb
    foreach d $dcdSet {
	set dcd $dir/nw_${name}${d}.dcd
	if {![checkExistence $dcd]} {set okay 0}
    }
}

if {!$okay} {exit}

# Analyze.
source $analysisScript
foreach s $sim {
    foreach {name moiety structPrefix dir dcdSet dcdFreq} $s {break}
    compute $name $moiety $structPrefix $dir $dcdSet $dcdFreq $outDir $probePdb $probeName
}

exit



