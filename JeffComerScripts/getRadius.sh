for typ in `cat dna_types.txt`; do
    grep "^${typ} " par_all27_na_nonbond.prm | awk '{print $1,$4}'
done
