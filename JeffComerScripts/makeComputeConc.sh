structPre=../dcd/nw_extra_layer_
vmdCmd="vmd -dispdev text -e concGridTraj.tcl -args"
stride=100
outDir=grid


for sys in pos neg nng ngn; do    
    if [ $sys == ngn ]; then
	dcdGlob="../dcd/nw_layer0*${sys}.dcd"
    else
	dcdGlob="../dcd/nw_layer{2,3a,4}*${sys}.dcd"
    fi

    struct=${structPre}${sys}
    name=dopc_${sys}
    dcdList=`eval echo $dcdGlob`
    
    echo "$vmdCmd $name POT $struct $outDir $stride $dcdList"
    echo "$vmdCmd $name CLA $struct $outDir $stride $dcdList"
done
