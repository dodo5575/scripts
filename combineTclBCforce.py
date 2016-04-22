#!/usr/bin/env python
# This script combines the tclBC force output from different patch.
# Usage: python combineTclBCforce.py input 
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
s="_combine.dat"
outputPrefix = inputPrefix[0]+s


## Calculate and combine
force = []
instantF = 0

for i in range(0,(len(array)-1)):

    if array[i][0] == array[i+1][0]:
        instantF = instantF + array[i][1]
    else:
        instantF = instantF + array[i][1]
        force.append([array[i][0],instantF])
        instantF = 0

instantF = instantF + array[i+1][1]
force.append([array[i+1][0],instantF])
instantF = 0


## Write output
output = open(outputPrefix,'w')

for i in range(0,len(force)):
    output.write('%f\t%f\n' % (force[i][0],force[i][1]))

output.close()
