
for f in lego_approach_*.z.dat; do
    val=$(awk 'BEGIN {sum=0; count=0}; $1 > -29 && $1 < -25 {sum+=$2; count++;}; END {print sum/count}' $f)

    shift=$(perl -e "print -($val)")
    echo $shift
    awk "{print \$1,\$2+$shift}" $f > shift_${f}
done
