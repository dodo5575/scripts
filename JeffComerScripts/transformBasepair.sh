cou0=566.443209

for diel in 80 92 114 124 142 176; do
    for sys in TA.100 TA.1000; do
	ion=pot
	cou=$(echo "$cou0/$diel" | bc -l 2> /dev/null)

	name=bp_diel${diel}_${sys}_${ion}
	f=transform_basepair/diel${diel}_scaled_distro_${sys}_${ion}.txt
	echo "gridTransformAddCoulomb $f $cou grids/$name.dx"

	name=bpx_diel${diel}_${sys}_${ion}
	f=transform_basepair/x_diel${diel}_scaled_distro_${sys}_${ion}.txt
	echo "gridTransformAddCoulomb $f $cou grids/$name.dx"
    done
done
