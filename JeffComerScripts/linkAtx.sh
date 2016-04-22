for f in pmf_diel*_at_chl.dx; do
    out=$(echo $f | sed 's/at/atx/')
    
    ln -s $f $out
done