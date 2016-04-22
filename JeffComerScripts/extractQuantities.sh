# Convert AA to nm.

for f in cyl_*.dat; do
    suf=${f#cyl_}
    
    awk '{print $1,$2,$3}' $f > num_${suf}; # in 1
    awk '{print $1,$4,$5}' $f > disp_${suf}; # in nm
    awk '{print $1,$6,$7}' $f > vel_${suf}; # in nm/ns
    awk '{print $1,$8,$9}' $f > conc_${suf}; # in M
    awk '{print $1,$10,$11}' $f > diffuse_${suf}; # in nm^2/ns
    awk '{print $1,$12,$13}' $f > curr_${suf}; # in pA
    awk '{print $1,$14,$15}' $f > j_${suf}; # in pA/nm^2
done

for f in curr_POT*.dat; do
    suf=${f#curr_POT_}

    paste curr_POT_${suf} curr_CLA_${suf} | awk '{print $1,$2+$5,sqrt($3*$3+$6*$6)}' > currBoth_${suf}
    awk '{sum+=$2;} END {print sum}' currBoth_${suf} > currSum_${suf}
done
