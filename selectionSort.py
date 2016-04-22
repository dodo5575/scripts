#!/usr/bin/env python
# This script sorts a (N x 2) array.
# Usage: python selectionSort.py input 
# Chen-Yu Li     cli56@illinois.edu
# 2014/3/5

## Standard module
import sys, os
import numpy as np
import re

## User's own module
sys.path.append('/home/cli56/scripts')
import ReadData

## Read input
print "Start loading data..."
array=ReadData.loadAscii(sys.argv[1])
print "Finished loading data!"

## set output file name
inputPrefix = re.split('\.', sys.argv[1])
inputPrefix_noType = ''
for i in range(len(inputPrefix)-1):
    inputPrefix_noType = inputPrefix_noType + inputPrefix[i]
    if i < (len(inputPrefix)-2):
        inputPrefix_noType = inputPrefix_noType + '.' 
    

s="_sorted.dat"
outputPrefix = inputPrefix_noType+s


## Compare and sort
print "Start sorting..."
for j in range(0,(len(array)-1)):

    imin = j

    for i in range(j+1,len(array)):

        if array[i][0] < array[imin][0]:

            imin = i

    if imin != j:

        tmp_value1 = array[imin][0]
        tmp_value2 = array[imin][1]

        array[imin][0] = array[j][0]
        array[imin][1] = array[j][1]

        array[j][0] = tmp_value1
        array[j][1] = tmp_value2
print "Finished sorting!"

## Write output
output = open(outputPrefix,'w')

for i in range(0,len(array)):
    output.write('%f\t%f\n' % (array[i][0],array[i][1]))

print "Finished writing output!"
output.close()
