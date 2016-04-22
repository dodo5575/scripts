set script pmeTraj.tcl
set name water
set moiety none
set structPrefix trap2.0_water
set dcd {output/trap2.0_water_1V1.dcd}
set dcdFreq 500
set outDir pme
set stride 1

source $script
compute $name $moiety $structPrefix $dcd $dcdFreq $outDir $stride

exit
