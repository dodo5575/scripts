# jcomer2@uiuc.edu
# Analysis:

set analysisScript currentThroughTraj1.tcl
set outDir current_through
#set outDir new

set outputPeriod 2000
set timestep 2e-5
set dcdFreq0 [expr 1e6*$outputPeriod*$timestep]
#set dcdFreq1 [expr 1e6*$outputPeriod*$timestep*0.5]

set set10 {0 1 2 3 4 5 6 7 8 9}
set set15 {0 1 2 3 4 5 6 7 8 9 10 11 12 13 14}
set set20 {0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19}
set set20a {0 1 2 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19}
set set20b {0 1 2 4 5 6 7 8 9 10 11 12 13 14 15 16 17}
set set4 {0 1 2 3}
set set4a {"" 0 1}
set set6a {1 2 3 4 5}
set set8 {0 1 2 3 4 5 6 7}
set set6 {0 1 2 3 4 5}
set set2 {0 1}

# Input data:
set sim {}
#lappend sim [list best_at_neg condor brown . $set6  $dcdFreq0 20]
#lappend sim [list best_at_pos condor brown . $set6  $dcdFreq0 20]
#lappend sim [list best_gc_neg condor brown . $set6  $dcdFreq0 20]
#lappend sim [list best_gc_pos condor brown . $set6  $dcdFreq0 20]
#lappend sim [list best_gcm_neg condor brown . $set6  $dcdFreq0 20]
#lappend sim [list best_gcm_pos condor brown . $set6  $dcdFreq0 20]
#lappend sim [list best_no_neg condor brown . $set6  $dcdFreq0 20]
#lappend sim [list best_no_pos condor brown . $set6  $dcdFreq0 20]

#lappend sim [list big_gc_neg condor brown . $set6  $dcdFreq0 20]
#lappend sim [list big_gc_pos condor brown . $set6  $dcdFreq0 20]

#lappend sim [list lowConc_at_neg condor brown . $set15  $dcdFreq0 20]
#lappend sim [list lowConc_at_pos condor brown . $set15 $dcdFreq0 20]
#lappend sim [list lowConc_gc_neg condor brown . $set15  $dcdFreq0 20]
#lappend sim [list lowConc_gc_pos condor brown . $set15  $dcdFreq0 20]

#lappend sim [list lowConc1_at_neg condor brown . $set20  $dcdFreq0 20]
#lappend sim [list lowConc1_at_pos condor brown . $set20a $dcdFreq0 20]
#lappend sim [list lowConc1_gc_neg condor brown . $set20  $dcdFreq0 20]
#lappend sim [list lowConc1_gc_pos condor brown . $set20  $dcdFreq0 20]
#lappend sim [list lowConc1_gc_neg condor brown . $set20  $dcdFreq0 20]
#lappend sim [list lowConc1_gc_pos condor brown . $set20  $dcdFreq0 20]

#lappend sim [list lowConc2_at_neg condor brown . $set20a  $dcdFreq0 20]
#lappend sim [list lowConc2_at_pos condor brown . $set20 $dcdFreq0 20]
#lappend sim [list lowConc2_gc_neg condor brown . $set20b  $dcdFreq0 20]
#lappend sim [list lowConc2_gc_pos condor brown . $set20  $dcdFreq0 20]
#lappend sim [list lowConc2_no_neg condor brown . $set20  $dcdFreq0 20]
#lappend sim [list lowConc2_no_pos condor brown . $set20  $dcdFreq0 20]

#lappend sim [list lowConc3_at_neg condor brown . $set20a  $dcdFreq0 20]
#lappend sim [list lowConc3_at_pos condor brown . $set20 $dcdFreq0 20]
#lappend sim [list lowConc3_gc_neg condor brown . $set20b  $dcdFreq0 20]
#lappend sim [list lowConc3_gc_pos condor brown . $set20  $dcdFreq0 20]
#lappend sim [list lowConc3_no_neg condor brown . $set20  $dcdFreq0 20]
#lappend sim [list lowConc3_no_pos condor brown . $set20  $dcdFreq0 20]

#lappend sim [list best1_at_neg condor brown . $set6  $dcdFreq0 20]
#lappend sim [list best1_at_pos condor brown . $set6  $dcdFreq0 20]
#lappend sim [list best1_gc_neg condor brown . $set6  $dcdFreq0 20]
#lappend sim [list best1_gc_pos condor brown . $set6  $dcdFreq0 20]
#lappend sim [list best1_gcm_neg condor brown . $set6  $dcdFreq0 20]
#lappend sim [list best1_gcm_pos condor brown . $set6  $dcdFreq0 20]
#lappend sim [list best1_no_neg condor brown . $set6  $dcdFreq0 20]
#lappend sim [list best1_no_pos condor brown . $set6  $dcdFreq0 20]

#lappend sim [list best2_at_neg condor brown . $set4  $dcdFreq0 20]
#lappend sim [list best2_at_pos condor brown . $set4  $dcdFreq0 20]
#lappend sim [list best2_gc_neg condor brown . $set4  $dcdFreq0 20]
#lappend sim [list best2_gc_pos condor brown . $set4  $dcdFreq0 20]
#lappend sim [list best2_gcm_neg condor brown . $set4  $dcdFreq0 20]
#lappend sim [list best2_gcm_pos condor brown . $set4  $dcdFreq0 20]

#lappend sim [list best3_at_neg condor brown . $set6a  $dcdFreq0 20]
#lappend sim [list best3_at_pos condor brown . $set6  $dcdFreq0 20]
#lappend sim [list best3_gc_neg condor brown . $set6  $dcdFreq0 20]
#lappend sim [list best3_gc_pos condor brown . $set6  $dcdFreq0 20]
#lappend sim [list best3_gcm_neg condor brown . $set6  $dcdFreq0 20]
#lappend sim [list best3_gcm_pos condor brown . $set6  $dcdFreq0 20]
#lappend sim [list best3_no_neg condor brown . $set6  $dcdFreq0 20]
#lappend sim [list best3_no_pos condor brown . $set6  $dcdFreq0 20]

lappend sim [list current_test test brown . {""} [expr 1e6*1*2e-5] 0]

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
    foreach {name moiety structPrefix dir dcdSet dcdFreq startFrame} $s {break}    

    checkExistence $structPrefix.xsc
    set dcdList {}
    foreach d $dcdSet {
	if {[string equal $moiety condor]} {
	    set dcd $dir/${name}.brown.$d.pdb
	} else {
	    set dcd $dir/${name}${d}.pdb
	}
	lappend dcdList $dcd
	if {![checkExistence $dcd]} {set okay 0}
    }
    lappend simDcdList $dcdList
}

if {!$okay} {exit}

# Analyze.
source $analysisScript
foreach s $sim dcdList $simDcdList {
    foreach {name moiety structPrefix dir dcdSet dcdFreq startFrame} $s {break}
    if {[string equal $moiety reg]} { set name ${name}0 }
    puts "Simulation: $s"
    compute $name $moiety $structPrefix $dcdList $dcdFreq $outDir $startFrame
}

exit
