#!/usr/bin/env python
# This script separate the data by quadrants and sort by density
# Usage: python getDNApeak.py input 
# Chen-Yu Li	cli56@illinois.edu
# 2014/9/3

import math, random, sys, os, re, matplotlib
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import cm
from numpy.random import randn


data = np.loadtxt(sys.argv[1])

## separate by axis
data_sep0 = np.array([])
data_sep1 = np.array([])
data_sep2 = np.array([])
data_sep3 = np.array([])

for i in range(len(data)):
    if data[i][0] >= 0 and data[i][1] >= 0:
        data_sep0 = np.append(data_sep0, data[i])
    if data[i][0] <  0 and data[i][1] >= 0:
        data_sep1 = np.append(data_sep1, data[i])
    if data[i][0] >= 0 and data[i][1] <  0:
        data_sep2 = np.append(data_sep2, data[i])
    if data[i][0] <  0 and data[i][1] <  0:
        data_sep3 = np.append(data_sep3, data[i])


data_sep0_reshape = data_sep0.reshape((data_sep0.size/5),5) 
data_sep1_reshape = data_sep1.reshape((data_sep1.size/5),5) 
data_sep2_reshape = data_sep2.reshape((data_sep2.size/5),5) 
data_sep3_reshape = data_sep3.reshape((data_sep3.size/5),5) 


## sort by density
data_sep0_reshape_sort = data_sep0_reshape[np.argsort(-data_sep0_reshape[:,2])]
data_sep1_reshape_sort = data_sep1_reshape[np.argsort(-data_sep1_reshape[:,2])]
data_sep2_reshape_sort = data_sep2_reshape[np.argsort(-data_sep2_reshape[:,2])]
data_sep3_reshape_sort = data_sep3_reshape[np.argsort(-data_sep3_reshape[:,2])]


(x0,y0) = data_sep0_reshape_sort.shape
(x1,y1) = data_sep1_reshape_sort.shape
(x2,y2) = data_sep2_reshape_sort.shape
(x3,y3) = data_sep3_reshape_sort.shape

## set output file name
inputPrefix = re.split('\.', sys.argv[1])
inputPrefix_noType = ''
for i in range(len(inputPrefix)-1):
    inputPrefix_noType = inputPrefix_noType + inputPrefix[i]
    if i < (len(inputPrefix)-2):
        inputPrefix_noType = inputPrefix_noType + '.' 

## write output

outname = inputPrefix_noType+"_getDNApeak.dat"
out = open(outname, "w")

out.write("#Top 5\n")
out.write("#        X         Y   Density\n")

for i in range(5):
    for j in range(y0):
        out.write("%10.3f" % (data_sep0_reshape_sort[i][j]))
    out.write("\n")
    
out.write("\n")
for i in range(5):
    for j in range(y1):
        out.write("%10.3f" % (data_sep1_reshape_sort[i][j]))
    out.write("\n")

out.write("\n")
for i in range(5):
    for j in range(y2):
        out.write("%10.3f" % (data_sep2_reshape_sort[i][j]))
    out.write("\n")

out.write("\n")
for i in range(5):
    for j in range(y3):
        out.write("%10.3f" % (data_sep3_reshape_sort[i][j]))
    out.write("\n")


out.write("\n")
out.write("#All\n")
out.write("\n")


for i in range(x0):
    for j in range(y0):
        out.write("%10.3f" % (data_sep0_reshape_sort[i][j]))
    out.write("\n")
    
out.write("\n")
for i in range(x1):
    for j in range(y1):
        out.write("%10.3f" % (data_sep1_reshape_sort[i][j]))
    out.write("\n")

out.write("\n")
for i in range(x2):
    for j in range(y2):
        out.write("%10.3f" % (data_sep2_reshape_sort[i][j]))
    out.write("\n")

out.write("\n")
for i in range(x3):
    for j in range(y3):
        out.write("%10.3f" % (data_sep3_reshape_sort[i][j]))
    out.write("\n")


out.close()


