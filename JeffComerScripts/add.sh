for f in ionK_*; do
    end=${f#*_}
    other=ionCl_${end}
    paste $f $other | awk {'print $1,$2+$4'} > ion_${end}
done
