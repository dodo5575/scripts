#!/usr/bin/env python
# This script generates colormap by averaging the dx file from a given axis.
# Usage: python plotDXaveColormap.py input dimension gridDist 
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

####################
## Load arguments ##
####################

infile = open(sys.argv[1])
theLines = infile.readlines()
line2 = theLines[1]
lastLine = theLines[-1] 
seg = line2[36:len(line2)]
dimArray = removeBlanksFromList(re.split('\s',seg))
for i in range(len(dimArray)):
    dimArray[i] = int(dimArray[i])

infile.close()
xDim = dimArray[0] 
yDim = dimArray[1] 
zDim = dimArray[2] 
dataNum = xDim*yDim*zDim 
dim = str(sys.argv[2])
gridDist = float(sys.argv[3])



########################
## Load input dx file ##
########################
data_raw = np.zeros((zDim,yDim,xDim))
if (dataNum % 3) == 0: 
    dxFile = np.loadtxt(sys.argv[1], delimiter=' ', skiprows=8)
    dxFile_flatten = dxFile.flatten() 
else:
    dxFile = np.genfromtxt(sys.argv[1], delimiter=' ', skip_header=8, skip_footer=1)
    lastLineList = removeBlanksFromList(re.split('\s',lastLine))
    dxFile_flatten = dxFile.flatten() 
    dxFile_flatten = np.append(dxFile_flatten,lastLineList)

#print dxFile


## set output file name
inputPrefix = re.split('\.', sys.argv[1])
inputPrefix_noType = ''
for i in range(len(inputPrefix)-1):
    inputPrefix_noType = inputPrefix_noType + inputPrefix[i]
    if i < (len(inputPrefix)-2):
        inputPrefix_noType = inputPrefix_noType + '.' 
    

s='_%s.pdf' % (dim)
outputName = inputPrefix_noType+s


################################################
## Transform input and save in to a 3-D array ##
################################################
tmp = 0
for i in range(xDim):
    for j in range(yDim):
        for k in range(zDim):

            data_raw[k][j][i] = dxFile_flatten[tmp]
            tmp += 1

#print data_raw
data = np.array([])
for k in range(zDim):
    if k == 0:
        data = data_raw[zDim-1].copy()
    else:
        data = np.append(data,data_raw[zDim-1-k],axis=0)

data = np.resize(data,(zDim,yDim,xDim))
#print data 


#################################
## Average along x,y or z axis ##
#################################
if dim == 'x': 
    data_ave = np.zeros((zDim,yDim))
    for k in range(zDim):
        for j in range(yDim):
    
            total = 0
            for i in range(xDim):
                total = total + data[k][j][i]
    
            total /= xDim
    
            data_ave[k][j] = total

elif dim == 'y':
    data_ave = np.zeros((zDim,xDim))
    for k in range(zDim):
        for i in range(xDim):
    
            total = 0
            for j in range(yDim):
                total = total + data[k][j][i]
    
            total /= yDim
    
            data_ave[k][i] = total

elif dim == 'z':
    data_ave = np.zeros((yDim,xDim))
    for j in range(yDim):
        for i in range(xDim):
    
            total = 0
            for k in range(zDim):
                total = total + data[k][j][i]
    
            total /= zDim
    
            data_ave[j][i] = total

#print data_ave

###########################
## Settings for plotting ##
###########################

#data_min = round(data_ave.min(),3)
#data_max = round(data_ave.max(),3)
data_min = trunc(data_ave.min(),2)
data_max = data_ave.max()
if data_max < 0:
    data_max = trunc((math.floor(data_max * 100) / 100.0),2)
else:
    data_max = trunc(data_ave.max(),2)


cbar_labels = np.linspace(data_min,data_max,num=5)  
for i in range(len(cbar_labels)):
    cbar_labels[i] = round(cbar_labels[i],2)  
cbar_labels_text = []
for i in range(len(cbar_labels)):
    #if i == (len(cbar_labels) - 1):
    if i == 0:
        cbar_labels_text.append(str(cbar_labels[i])+' (V)')  
    else:       
        cbar_labels_text.append(str(cbar_labels[i]))  


#print data_min
#print data_max


###########################
## Plotting the colormap ##
###########################

# Make plot with vertical (default) colorbar
fig, ax = plt.subplots()

if dim == 'x':

    cax = ax.imshow(data_ave, extent=[0,yDim*gridDist,0,zDim*gridDist], interpolation='nearest', cmap=cm.jet)
    #ax.set_title('Electrostatic Potential', fontsize=32)
    xlabel = ax.set_xlabel(r'Y ($\AA$)', fontsize=28)
    ylabel = ax.set_ylabel(r'Z ($\AA$)', fontsize=28)
    ax.tick_params(labelsize=24) 
    xTickInt = int(int(yDim*gridDist/10)/4) * 10
    yTickInt = int(int(zDim*gridDist/10)/4) * 10

elif dim == 'y':

    cax = ax.imshow(data_ave, extent=[0,xDim*gridDist,0,zDim*gridDist], interpolation='nearest', cmap=cm.jet)
    #ax.set_title('Electrostatic Potential', fontsize=32)
    xlabel = ax.set_xlabel(r'X ($\AA$)', fontsize=28)
    ylabel = ax.set_ylabel(r'Z ($\AA$)', fontsize=28)
    ax.tick_params(labelsize=24) 
    xTickInt = int(int(xDim*gridDist/10)/4) * 10
    yTickInt = int(int(zDim*gridDist/10)/4) * 10

elif dim == 'z':

    cax = ax.imshow(data_ave, extent=[0,xDim*gridDist,0,yDim*gridDist], interpolation='nearest', cmap=cm.jet)
    #ax.set_title('Electrostatic Potential', fontsize=32)
    xlabel = ax.set_xlabel(r'X ($\AA$)', fontsize=28)
    ylabel = ax.set_ylabel(r'Y ($\AA$)', fontsize=28)
    ax.tick_params(labelsize=24) 
    xTickInt = int(int(xDim*gridDist/10)/4) * 10
    yTickInt = int(int(yDim*gridDist/10)/4) * 10


# set tick size
ax.tick_params('both', length=10, width=1.5, which='major')
ax.tick_params('both', length=5, width=1, which='minor')
# set axis width
for axis in ['top','bottom','left','right']:
  ax.spines[axis].set_linewidth(1.5)

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

# Add colorbar, make sure to specify tick locations to match desired ticklabels
#cbar = fig.colorbar(cax, ticks=[data_min, data_max])

# create an axes on the right side of ax. The width of cax will be 5%
# of ax and the padding between cax and ax will be fixed at 0.05 inch.
divider = make_axes_locatable(ax)
caax = divider.append_axes("right", size="5%", pad=0.05)

cbar = fig.colorbar(cax, cax=caax)
#cbar.set_label('mV', fontsize=24)# vertically oriented colorbar
cbar.set_ticks(cbar_labels)
cbar.ax.set_yticklabels(cbar_labels_text,fontsize=24)# vertically oriented colorbar

# set tick size
cbar.ax.tick_params('both', length=5, width=1.5, which='major')
# set axis width
for axis in ['top','bottom','left','right']:
  cbar.ax.spines[axis].set_linewidth(1.5)
cbar.ax.invert_yaxis()

#cbar.ax.set_yticklabels([data_min, data_max], fontsize=24)# vertically oriented colorbar

#plt.savefig('plot.pdf', bbox_extra_artists=[xlabel], bbox_inches='tight')

plt.savefig(outputName, bbox_inches='tight')
plt.show()

