#!/usr/bin/env python
# This script calculates mean area from equilibration dimension in the last 400ns and find the closest restart file
# Usage: python findCloseRestart.py EQdim dimension xscPrefix numOfXsc outName
# Chen-Yu Li     cli56@illinois.edu
# 2014/6/6


#import stats
import math, random, sys, os, re, matplotlib 
import numpy as np
import matplotlib.pyplot as plt
#from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import cm
from numpy.random import randn

## User's own module
sys.path.append('/home/cli56/scripts')
import ReadData
import selSort


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
## Load parameters  ####
########################

EQdim = sys.argv[1]
dim = str(sys.argv[2])
xscPrefix = sys.argv[3]
numOfXsc = int(sys.argv[4])
outName = sys.argv[5] 


########################
## Load input file  ####
########################
array_raw = np.loadtxt(EQdim, delimiter=' ')
N_raw = len(array_raw) 
## load only the last 400 ns 
array = array_raw[(N_raw - 166667):N_raw, : ] 
N = len(array)


## set output file name
s='_last400ns_findCloseRestart.dat'
outputName = outName+s


## set xsc list
xscList = []
for i in range(numOfXsc):
    s = ".%d.restart.xsc" % (i+1)
    xscList.append(xscPrefix+s)

#print xscList
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


## get average dimension
meanX = np.mean(array[:,1]) 
meanY = np.mean(array[:,2]) 
meanZ = np.mean(array[:,3]) 
meanArea = np.mean(area)
meanXYratio = meanX/meanY


## load data from all xsc files
xsc_array = np.array([])
for i in range(numOfXsc):
    xsc_array_n = np.loadtxt(xscList[i],skiprows=2)
    xsc_array = np.append(xsc_array,xsc_array_n)

xsc_array = np.resize(xsc_array,(numOfXsc,19))
#print xsc_array

## calculate area and XY ratio in the xsc files
xscArea = []
dArea = np.array([])
dXYratio = np.array([])
for i in range(numOfXsc):
    xscArea_n = xsc_array[i][1] * xsc_array[i][5]
    xscArea.append(xscArea_n)
    dArea = np.append(dArea,[i,abs(meanArea - xscArea_n)])
    #print [i,meanXYratio, xsc_array[i][1], xsc_array[i][5]]
    dXYratio = np.append(dXYratio,[i,abs(meanXYratio - xsc_array[i][1] / xsc_array[i][5])])


dArea = np.resize(dArea,(numOfXsc,2))
dXYratio = np.resize(dXYratio,(numOfXsc,2))
dArea_sort = selSort.selSort(dArea,1)
dXYratio_sort = selSort.selSort(dXYratio,1)
#print dArea_sort
#print dXYratio_sort


###################
## write output####
###################
output = open(outputName, "w")
output.write("EQ Dim mean x y z area XYratio = %.3f\t%.3f\t%.3f\t%.3f\t%.3f\n\n" % (meanX,meanY,meanZ,meanArea,meanXYratio))

output.write("name\tdArea\tarea\tx\ty\tz\n")
for tuple in dArea_sort:
    output.write("%s\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n" % (xscList[int(tuple[0])], tuple[1], xscArea[int(tuple[0])], xsc_array[int(tuple[0])][1], xsc_array[int(tuple[0])][5], xsc_array[int(tuple[0])][9]))
output.write("\n")

output.write("name\tdXYratio\tXYratio\tx\ty\tz\n")
for tuple in dXYratio_sort:
    output.write("%s\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n" % (xscList[int(tuple[0])], tuple[1], (xsc_array[int(tuple[0])][1] / xsc_array[int(tuple[0])][5]), xsc_array[int(tuple[0])][1], xsc_array[int(tuple[0])][5], xsc_array[int(tuple[0])][9]))
output.write("\n")


output.close()
