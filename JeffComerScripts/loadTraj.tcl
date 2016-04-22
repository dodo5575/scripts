set name dcd/nw_bam_specific_1V
set dcdList {0.0 0.1 1.0 1.1 2.0 2.1 3.0 3.1 4.0 4.1 5.0 5.1 6.0 6.1 7.0 7.1 8.0 8.1 9.0 9.1 10.0 10.1 11.0}

#set name dcd/nw_bam_specific_1.5V 
#set dcdList {0.0 0.1 1.0 1.1 2.0 3.0 3.1 3.2 4.0 4.1 5.0 5.1 5.2 6.0 6.1 7a.0 8.0 8.1 8.2}

set name dcd/nw_bam_specific_2V
set dcdList {0.0 1.0 1.1 2.0 2.1 3.0 3.1}

foreach dcd $dcdList {
    mol addfile ${name}${dcd}.dcd step 20 waitfor all
}
