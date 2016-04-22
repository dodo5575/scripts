set sys 3
mol load psf nw_bam_spec_pull${sys}.psf

set dcdDir /scratch/tbgl/jcomer/protein-dna/pull/dcd
set dcdList [lsort [glob $dcdDir/nw_specbam_pull_smd*_sys${sys}.dcd]]

foreach dcd $dcdList {
    mol addfile ${dcd} step 20 waitfor all
}

source representComplexCartoon.tcl
