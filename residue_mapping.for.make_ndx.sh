# convert the residues in the make_ndx files using a residue mapping file
# Usage: ./residue_mapping.for.make_ndx.sh map input output
# Author: Chen-Yu Li <cli56@illinois.edu> 
# 2015/11/21


map=$1
in=$2
out=$3

rm -f $out

declare -A Array

numberoflinestoskip=1
mapfile -s $numberoflinestoskip -t myArray < $map


regex='([0-9]*)\s*([0-9]*)'
for i in ${!myArray[@]}
do

    [[ ${myArray[$i]} =~ $regex ]]
    Array[${BASH_REMATCH[1]}]=${BASH_REMATCH[2]}

done


regex='ri\s*([0-9]*)\s*([0-9]*)'

while read -ra line
do

    if [[ ${line[@]} =~ $regex ]]
    then

        # convert to zero-based, as in vmd
        pre_r1=`echo "${BASH_REMATCH[1]} - 1" | bc -l`

        if [ ${Array[$pre_r1]} ]
        then
            # convert back to one-based for gromacs
            r1=`echo "${Array[$pre_r1]} + 1" | bc -l`
        else
            r1=${BASH_REMATCH[1]}
        fi

        if [ ${BASH_REMATCH[2]} ]
        then
            pre_r2=`echo "${BASH_REMATCH[2]} - 1" | bc -l`

            if [ ${Array[$pre_r2]} ]
            then
                r2=`echo "${Array[$pre_r2]} + 1" | bc -l`
            else
                r2=${BASH_REMATCH[2]}
            fi

            final=`echo ${line[@]} | sed "s/${BASH_REMATCH[1]}/$r1/" | sed "s/${BASH_REMATCH[2]}/$r2/"`

        else
            final=`echo ${line[@]} | sed "s/${BASH_REMATCH[1]}/$r1/"`
        fi


        echo $final >> $out 

    else
        echo "${line[@]}" >> $out

    fi   


done < $in 

echo "quit" >> $out




