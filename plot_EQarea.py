#!/usr/bin/env python
# This script plots area trace of equilibration.
# Usage: python plot_EQarea.py input dimension xstFrequency timestep
# Chen-Yu Li     cli56@illinois.edu
# 2014/5/8


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
array = np.loadtxt(sys.argv[1], delimiter=' ')
dim = str(sys.argv[2])

# interval (ns)
interval = float(sys.argv[3]) * float(sys.argv[4]) * 10**(-6)


## set output file name
inputPrefix = re.split('\.', sys.argv[1])
inputPrefix_noType = ''
for i in range(len(inputPrefix)-1):
    inputPrefix_noType = inputPrefix_noType + inputPrefix[i]
    if i < (len(inputPrefix)-2):
        inputPrefix_noType = inputPrefix_noType + '.'
        

s='_area.pdf'
outputName = inputPrefix_noType+s



########################################
## Calculate area along x,y or z axis ##
########################################
time = []
for i in range(len(array)):
    time.append(interval * i)


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


###########################
## Settings for plotting ##
###########################

#data_min = round(data_ave.min(),3)
#data_max = round(data_ave.max(),3)
#data_min = trunc(data_ave.min(),2)
#data_max = data_ave.max()
#if data_max < 0:
#    data_max = trunc((math.floor(data_max * 100) / 100.0),2)
#else:
#    data_max = trunc(data_ave.max(),2)
#
#
#cbar_labels = np.linspace(data_min,data_max,num=5)  
#for i in range(len(cbar_labels)):
#    cbar_labels[i] = round(cbar_labels[i],2)  
#cbar_labels_text = []
#for i in range(len(cbar_labels)):
#    if i == (len(cbar_labels) - 1):
#        cbar_labels_text.append(str(cbar_labels[i])+' (V)')  
#    else:       
#        cbar_labels_text.append(str(cbar_labels[i]))  


#print data_min
#print data_max


###########################
## Plotting ###############
###########################

fig, ax = plt.subplots()
plt.plot(time,area)
plt.xlabel('Time (ns)', fontsize=28)
plt.ylabel(r'Area ($\AA$)', fontsize=28)
ax.tick_params(labelsize=24)
ax.yaxis.grid()

xTickInt = 50
yTickInt = 100


majorXLocator   = matplotlib.ticker.MultipleLocator(xTickInt*2)
minorXLocator   = matplotlib.ticker.MultipleLocator(xTickInt)
majorXFormatter = matplotlib.ticker.FormatStrFormatter('%d')
ax.xaxis.set_major_locator(majorXLocator)
ax.xaxis.set_minor_locator(minorXLocator)
ax.xaxis.set_major_formatter(majorXFormatter)

majorYLocator   = matplotlib.ticker.MultipleLocator(yTickInt*2)
minorYLocator   = matplotlib.ticker.MultipleLocator(yTickInt)
majorYFormatter = matplotlib.ticker.FormatStrFormatter('%d')
ax.yaxis.set_major_locator(majorYLocator)
ax.yaxis.set_minor_locator(minorYLocator)
ax.yaxis.set_major_formatter(majorYFormatter)



plt.savefig(outputName, bbox_inches='tight')
plt.show()

