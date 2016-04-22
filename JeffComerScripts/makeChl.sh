for f in *.trans; do
    out=$( echo $f | sed 's/_pot/_chl/' )
    sed 's/_pot/_chl/' < $f > $out
done
