function findBadNodes() 
{
    out=ref_mean_$1.dat
    echo -n "" > $out
    for f in *$1*.log; do
	ind=`grep "INDEX: " $f`
	ref=`grep "REF: " $f`
	mean=`grep "MEAN: " $f`

	echo "$ind $ref $mean" >> $out
    done

# Filter for the bad nodes.
    echo "Finding bad nodes..."
    awk '{print $2,($4-$8)*($4-$8) + ($5-$9)*($5-$9) + ($6-$10)*($6-$10)}' $out | awk '($2>0.2) {print $1,$2}' > bad_nodes_$1.dat
}

findBadNodes at_pot
findBadNodes at_1M_pot
findBadNodes at_chl
findBadNodes gc_pot
