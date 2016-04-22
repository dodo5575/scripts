for f in pmf_ready*dx; do
    base=${f%.*}
    gridProfile $f x 0.5 0.5 $base.x.dat
done
