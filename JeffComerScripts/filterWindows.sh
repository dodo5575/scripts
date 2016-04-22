for f in window_*.dat; do
    base=${f%.*}
    ./checkWindowIndex.tcl $f $base.win
done
