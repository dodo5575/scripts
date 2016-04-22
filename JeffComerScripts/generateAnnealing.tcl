# Make an NAMD files for an annealing schedule.
# Schedule is {{temp0 dur0} {temp1 dur1}...}
set schedule {{0 500} {7000 20000} {5000 20000} {2000 100000} {300 50000}}
set structName silica_undulate_box
set lx "115 0 0"
set ly "0 115 0"
set lz "0 0 160"
set pmeGridSizeX 96
set pmeGridSizeY 96
set pmeGridSizeZ 128
set gridforceVFile exclude2-3-3_pore22.dx
set scale 50.0
# Output
set outName undulate22_anneal

set tab "\t\t\t"
set n 0
foreach s $schedule {
    set out [open ${outName}${n}.namd w]

    set temp [lindex $s 0]
    set steps [lindex $s 1]

    puts $out "\# NAMD config for annealing schedule $schedule"
    puts $out "numSteps${tab}${steps}"
    puts $out "structure${tab}${structName}.psf"
    puts $out "coordinates${tab}${structName}.pdb"
    puts $out "set output output/${outName}$n"
    puts $out "set last output/${outName}[expr $n-1]"

    puts $out ""
    puts $out "outputName${tab}\$output"
    if {$n == 0} {
	puts $out "temperature${tab}${temp}"
	puts $out "cellBasisVector1${tab}${lx}"
	puts $out "cellBasisVector2${tab}${ly}"
	puts $out "cellBasisVector3${tab}${lz}"
	
    } else {
	puts $out "binCoordinates${tab}\$last.coor"
	puts $out "binVelocities${tab}\$last.vel"
	puts $out "extendedSystem${tab}\$last.xsc"
    }
    
    if {$temp != 0} {
	puts $out ""
	puts $out "\# temperature control"
	puts $out "langevin${tab}on"
	puts $out "langevinTemp${tab}${temp}"
	puts $out "langevinDamping${tab}5"
	puts $out "langevinHydrogen${tab}off"
    } else {
	puts $out "minimization${tab}on"
    }

    puts $out ""
    puts $out "\# parameters"
    puts $out "paraTypeCharmm${tab}on"
    puts $out "parameters${tab}SiOtab.par"
    puts $out "tabulatedEnergies${tab}on"
    puts $out "tabulatedEnergiesFile${tab}bkstab.dat"
    puts $out "tableInterpType${tab}cubic"

    puts $out ""
    puts $out "exclude${tab}scaled1-4"
    puts $out "1-4scaling${tab}0.83333333"

    puts $out ""
    puts $out "switching${tab}on"
    puts $out "switchDist${tab}5.4"
    puts $out "cutoff${tab}5.5"
    puts $out "pairListDist${tab}8"
    puts $out "margin${tab}3"


    puts $out ""
    puts $out "\# integration"
    puts $out "timestep${tab}1"

    puts $out ""
    puts $out "\# output"
    puts $out "binaryOutput${tab}on"
    puts $out "binaryRestart${tab}on"
    puts $out "wrapAll${tab}on"
    puts $out "wrapNearest${tab}on"
    
    puts $out ""
    puts $out "outputEnergies${tab}1000"
    puts $out "outputTiming${tab}1000"
    puts $out "xstFreq${tab}5000"
    puts $out "dcdFreq${tab}5000"
    puts $out "restartFreq${tab}5000"

    puts $out ""
    puts $out "\# electrostatics"
    puts $out "pme${tab}on"
    puts $out "pmeGridSizeX${tab}${pmeGridSizeX}"
    puts $out "pmeGridSizeY${tab}${pmeGridSizeY}"
    puts $out "pmeGridSizeZ${tab}${pmeGridSizeZ}"
    puts $out "zeroMomentum${tab}on"
    
    puts $out ""
    puts $out "\# external forces"
    puts $out "gridforce${tab}on"
    puts $out "gridforceFile${tab}mark_all.pdb"
    puts $out "gridforceCol${tab}B"
    puts $out "gridforceQCol${tab}O"
    puts $out "gridforceVFile${tab}${gridforceVFile}"
    puts $out "gridforceScale${tab}${scale} $scale $scale"
    puts $out "gridforceCont1${tab}on"
    puts $out "gridforceCont2${tab}on"

    close $out
    incr n
}

set out [open $outName.csh w]
set n 0
foreach s $schedule {
    puts $out "namd2tab ${outName}${n}.namd >! ${outName}${n}.log"
    incr n
}
close $out



