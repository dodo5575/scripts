#!/bin/bash

inPut=$1

outPut="${inPut:0:$(expr ${#inPut} - 4)}_10H.psf"

sed 's/1\.0080/10\.080/g' $inPut > $outPut
sed -i 's/1\.0079/10\.079/g' $outPut


