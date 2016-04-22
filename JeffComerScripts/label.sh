for name in 100cM 140cM ""; do
    if [ -n "$name" ]; then
	f=$name
    else
	f=10cM
    fi


    sed s/\+// over_triplet${name}_unique.txt > tmp.dat

    awk '{printf "@    xaxis  tick major %s, %s\n@    xaxis  ticklabel %s, \"%s\"\n", NR-1, NR-1, NR-1, toupper($1)}' tmp.dat > label_triplet${f}.dat
    awk '{print NR-1,$2,$3}' tmp.dat > graph_triplet${f}.dat

    rm -f triplet_labeled_${f}.agr
    rm -f triplet_graph_${f}.agr

    tclsh templateGraph.tcl vertical_ticks.agr graph_triplet${f}.dat triplet_graph_${f}.agr
    tclsh labelGraph.tcl triplet_graph_${f}.agr label_triplet${f}.dat triplet_labeled_${f}.agr
done
