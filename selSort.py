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


## Function 
def selSort(inArray, col):
    array = inArray.copy()

    for j in range(0,(len(array)-1)):
    
        imin = j
    
        for i in range(j+1,len(array)):
    
            if array[i][col] < array[imin][col]:
    
                imin = i

        if imin != j:
    
            tmp_value1 = array[imin][0]
            tmp_value2 = array[imin][1]
    
            array[imin][0] = array[j][0]
            array[imin][1] = array[j][1]
    
            array[j][0] = tmp_value1
            array[j][1] = tmp_value2

    return array


