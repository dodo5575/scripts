# jcomer2@uiuc.edu
# Analysis:

#set analysisScript currentThroughTraj.tcl
#set outDir current

#set analysisScript permeationTraj.tcl
#set outDir perm

set analysisScript initialBasepairsTraj.tcl
set outDir basepair

#set analysisScript regionNumberTraj.tcl
#set outDir number

#set analysisScript positionTraj.tcl
#set analysisScript positionTrajProt.tcl
#set outDir pos
set stride 0

# Input data:
set sim {}
lappend sim {neutra_at dna ../nw_NeuA-pAT-pore_E ../output/nw_simulationAT {-1 -2} 5000}
lappend sim {neutra_gc dna ../nw_NeuA-pCG-pore_E ../output/nw_simulationCG {-1 -2} 5000}

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
    foreach {name moiety structPrefix dcdPrefix dcdSet dcdFreq} $s {break}    

    checkExistence $structPrefix.psf
    checkExistence $structPrefix.pdb
    checkExistence $structPrefix.xsc
    set dcdList {}
    foreach d $dcdSet {
	set dcd ${dcdPrefix}${d}.dcd
	lappend dcdList $dcd
	if {![checkExistence $dcd]} {set okay 0}
    }
    lappend simDcdList $dcdList
}

if {!$okay} {exit}

# Analyze.
source $analysisScript
foreach s $sim dcdList $simDcdList {
    foreach {name moiety structPrefix dir dcdSet dcdFreq} $s {break}
    puts "Simulation: $s"
    compute $name $moiety $structPrefix $dcdList $dcdFreq $outDir $stride
}

exit
