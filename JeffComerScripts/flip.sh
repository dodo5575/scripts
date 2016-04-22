#for f in *null_pos*_z-*dat; do
#    out=${f%_pos*}_neg_z${f##*_z-}
#    echo $f $out
#    awk '{print $1,-1.0*$2}' $f > $out
#done

#for f in *null_pos*_z1*dat; do
#    out=${f%_pos*}_neg_z-1${f##*_z1}
#    echo $f $out
#    awk '{print $1,-1.0*$2}' $f > $out
#done

for f in *null_pos*.dat; do
    out=${f%_pos*}_neg${f##*_pos}
    awk '{print $1,-1.0*$2}' $f > $out
done
