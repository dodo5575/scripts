for int in 1 2; do
    for dt in 2.5 5 10 20 40; do
	out=force_cylinder/cyl_force${int}_timestep${dt}.dat
	echo $out

	echo -n "" > $out
	for f in force${int}_timestep${dt}.brown.*.force; do
	    awk '($4 > -12 && $4 < 12) {print sqrt($2*$2+$3*$3),($2*$5 + $3*$6)/sqrt($2*$2+$3*$3)}' $f >>$out
	done
    done
done
