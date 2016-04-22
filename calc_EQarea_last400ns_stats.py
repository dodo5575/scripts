#!/usr/bin/env python
# This script calculates statistical properties from equilibration dimension in the last 400ns.
# Usage: python calc_EQarea_stats.py input dimension
# Chen-Yu Li     cli56@illinois.edu
# 2014/5/20


import stats
import math, random, sys, os, re, matplotlib 
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import cm
from numpy.random import randn


## Functions --------------------------------------------------------

def trunc(f, n):
    #'''Truncates/pads a float f to n decimal places without rounding'''
    slen = len('%.*f' % (n, f))
    return float(str(f)[:slen])


def removeBlanksFromList(theList):
    newList=[]
    for counter in range(len(theList)):
        if theList[counter]!='':
            newList.append(theList[counter])
    return newList


# ------------------------------------------------------------------------


########################
## Load input file ####
########################
array_raw = np.loadtxt(sys.argv[1], delimiter=' ')
N_raw = len(array_raw) 
array = array_raw[(N_raw - 166667):N_raw, : ] 
N = len(array)
dim = str(sys.argv[2])


## set output file name
inputPrefix = re.split('\.', sys.argv[1])
inputPrefix_noType = ''
for i in range(len(inputPrefix)-1):
    inputPrefix_noType = inputPrefix_noType + inputPrefix[i]
    if i < (len(inputPrefix)-2):
        inputPrefix_noType = inputPrefix_noType + '.'
        

s='_last400ns_stats.dat'
outputName = inputPrefix_noType+s



########################################
## Calculate area along x,y or z axis ##
########################################

if dim == 'x': 
    area = []
    for i in range(len(array)):
        area.append(array[i][2] * array[i][3])

elif dim == 'y':
    area = []
    for i in range(len(array)):
        area.append(array[i][1] * array[i][3])

elif dim == 'z':
    area = []
    for i in range(len(array)):
        area.append(array[i][1] * array[i][2])

result = stats.stats(area)

output = open(outputName, "w")
output.write("N = %d\n" % N)
output.write("mean = %f\n" % result[0])
output.write("stdev = %f\n" % result[1])
output.write("autocorrelation time = %f\n" % result[2])
output.write("error = %f\n" % result[3])
output.write("mean x y z = %f %f %f" % (np.mean(array[:,1]),np.mean(array[:,2]),np.mean(array[:,3])))

output.close()
