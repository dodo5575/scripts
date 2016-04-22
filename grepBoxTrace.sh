#!/bin/bash
#This scripts extracts the simulation box size information form xst file.
#Usage: ./grepBoxTrace.sh name xstPrefix startTime numOfxst
#By Chen-Yu Li	cli56@illinois.edu
#2014/2/2

for i in $(seq 1 ${4})
do

time=`echo "(${i}-1)*9.6+${3}" | bc -l`
input=${2}.${i}.xst
grep "0" $input | awk "{print \$1*0.000002+${time}, \$2}" >> ${1}_X.dat
grep "0" $input | awk "{print \$1*0.000002+${time}, \$6}" >> ${1}_Y.dat
grep "0" $input | awk "{print \$1*0.000002+${time}, \$10}" >> ${1}_Z.dat

done

