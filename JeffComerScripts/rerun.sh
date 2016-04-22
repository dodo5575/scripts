for config in nuc1_cyt_{pot,chl}_pos{0..340}.namd; do
    log=${config%.*}.log
    #echo $log

    if [ ! -e $log ]; then
	echo $config
    fi
done
