i=$1
j=$(echo $i | sed 's/\.pdb$/_NOATOMNUM.PDB/')
if [ -f $j ]; then
    echo "$j already exists!"
    exit -1
fi
awk 'BEGIN {OFS="\t"} /^ATOM/ {print $1, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12}' $i > $j
echo "Wrote file $j"