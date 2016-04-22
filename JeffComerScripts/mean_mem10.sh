t0=8
for file in dna_dna_mem10_1V_*.dat; do
    base=${file%.*}
    num=${base##*_}

    echo -n "$num "
    tclsh meanValueCut.tcl $t0 $file
    echo ""
done
