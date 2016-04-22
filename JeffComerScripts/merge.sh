for x in at_neg at_pos gc_neg gc_pos no_neg no_pos; do
    cat thru28_curr_mole_$x.dat thru28_curr_100mMa_$x.dat > thru28_curr_mergea_$x.dat
done

