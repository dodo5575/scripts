#!/usr/bin/env python
# This script invert the number in each column in a data set.
# Usage: python invert.py input 
# Chen-Yu Li     cli56@illinois.edu
# 2014/3/4

## Standard module
import sys, os
import numpy as np
import re

## User's own module
sys.path.append('/home/cli56/scripts')
import ReadData

## Read input
array=ReadData.loadAscii(sys.argv[1])


## set output file name
inputPrefix = re.split('\.', sys.argv[1])
s="_invert.dat"
outputPrefix = inputPrefix[0]+s


## Calculate and combine

for i in range(len(array)):
    for j in range(len(array[i])):
        if array[i][j] != 0:
            array[i][j] = 1.0 / array[i][j]


## Write output
output = open(outputPrefix,'w')

for i in range(0,len(array)):
    for j in range(len(array[i])):
        output.write('%e\t' % (array[i][j]))
        if j == (len(array[i]) - 1):
           output.write('\n') 


output.close()
