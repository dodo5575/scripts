function doSystem()
{
    local sys=$1
    local ion=$2

    out=diffNode_${sys}_${ion}.dat
    echo -n "" > $out
    for f in *${sys}_${ion}*.log; do

	index=`grep "INDEX: " $f | awk '{print $2}'`
	meanPos=`grep "MEAN: " $f | awk '{print $2,$3,$4}'`    
	diffuse=`grep "DIFFUSE: " $f | awk '{print $2,$3,$4}'`
	echo "$index $meanPos $diffuse" >> $out
    done

    eval "vmdtext -e makeDiffusionPdbFromNodes.tcl -args $out ../../run_${sys} ../diffusion/results_${sys}_${ion}"
}

doSystem at pot
doSystem at chl
#doSystem at_1M pot
#doSystem gc pot

