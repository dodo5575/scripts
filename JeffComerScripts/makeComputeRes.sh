structPre=../dcd/nw_pore40+dopc_
vmdCmd="vmd -dispdev text -e residueTrajectoryExtractor.tcl -args"
dcdGlob='../dcd/nw_eq{5,6}_dopc_${sys}.dcd ../dcd/nw_scram{0,1}_dopc_${sys}.dcd'
outDir=follow_res

for sys in pos neg neu; do
    for region in inTop inBot outTop outBot; do
	struct=${structPre}${sys}1
	name=dopc_${sys}
	dcdList=`eval echo "$dcdGlob"`
	
	echo "$vmdCmd $name $region $struct $outDir $dcdList"
    done
done
