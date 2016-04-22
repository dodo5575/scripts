for ion in pot chl; do
    for bp in ade thy gua; do
	name=pmf_pmfHard_${bp}_${ion}.dx
	mask=pmf_embed_${bp}_${ion}.dx.mask
	awk "{print \"$name\",\"$mask\",\$2,\$3,\$4,\$5,\$6,\$7,\$8,\$9,\$10,\$11,\$12,\$13}" trans_${bp}_${ion}.txt > transp_${bp}_${ion}.txt
    done
done
