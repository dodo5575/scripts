for sim in {0..19}; do
    for field in pos neg; do
	for sys in gcm gcm_1M; do
	    out=eq0_${sys}_${field}${sim}.namd 
	    echo $out

	    echo -n "" > $out
	    echo "set sys $sys" >> $out
	    echo "set field $field" >> $out
	    echo "set sim $sim" >> $out
	    echo "" >> $out
	    echo "source ../eq0.namd" >> $out
	done
    done
done
