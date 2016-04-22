# jcomer2@uiuc.edu
# Analysis:

# Parameters:
set analysisScript permeationTraj.tcl
set outDir permeation
set name trap_coated_2V
set moiety DNA
set structPrefix nw_coated_trap_charmm
set dcdPrefix dcd/nw_trap_coated_2V
set dcdId {0 1 2}
set dcdFreq 5000

# Make the list of dcd files.
set dcdList {}
foreach d $dcdId {
    lappend dcdList ${dcdPrefix}${d}.dcd
}

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
if {![checkExistence $structPrefix.psf]} {set okay 0}
if {![checkExistence $structPrefix.pdb]} {set okay 0}
foreach dcd $dcdList {
    if {![checkExistence $dcd]} {set okay 0}
}
if {!$okay} {exit}

# Analyze.
source $analysisScript
compute $name $moiety $structPrefix $dcdList $dcdFreq $outDir
exit
