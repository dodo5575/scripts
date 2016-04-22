for num in 4 30 40 50 60; do
    # A formula for estimating the required sim length.
    steps=$(perl -e "print int(2.0e10/($num*$num + 16))")

    sed -e "s/STEPS/$steps/" -e "s/NUM/$num/" < template_small_num.brown > small24_num${num}.brown

    sed -e "s/STEPS/$steps/" -e "s/NUM/$num/" -e 's/phantom24_40-40-72_soft.dx/bulk_40-40-72.dx/' < template_small_num.brown > bulk_num${num}.brown
done
