#!/usr/bin/env python
# This script generates colormap by averaging the dx file from a given axis.
# Usage: python plotStreamline.py input 


import math, random, sys, os, re, matplotlib
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
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



## Plotting Parameters ---------------------------------------------

titleSize = 32
labelSize = 28
axisSize = 24
xTickInt = 2.5
yTickInt = 2.5
cbarInt = 600
threshold = 100


####################
## Load arguments ##
####################

infile = open(sys.argv[1])
theLines = infile.readlines()

#for i, line in enumerate(infile):
#    if i == 3:
#        dimArray = removeBlanksFromList(re.split('\s',line[14:len(line)]))
#        print dimArray
#    elif i > 3:
#        break
#
#for i in range(len(dimArray)):
#    dimArray[i] = int(dimArray[i])

infile.close()

#xDim = dimArray[0]
#yDim = dimArray[1]
#zDim = dimArray[2]
#dataNum = xDim*yDim*zDim
dim = 'y'


## set output file name
inputPrefix = re.split('\.', sys.argv[1])
inputPrefix_noType = ''
for i in range(len(inputPrefix)-1):
    inputPrefix_noType = inputPrefix_noType + inputPrefix[i]
    if i < (len(inputPrefix)-2):
        inputPrefix_noType = inputPrefix_noType + '.'


s='_%s.pdf' % (dim)
outputName = inputPrefix_noType+s


########################
## Load input file #####
########################


regexp = r"([^\s]+)\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)"
output = np.fromregex(sys.argv[1], regexp,
	        [('x', np.float), ('z', np.float), ('dens', np.float), ('fx', np.float), ('fz', np.float) ])


X = np.reshape(output['x'],    (26,26))
Y = np.reshape(output['z'],    (26,26))
D = np.reshape(output['dens'], (26,26))
U = np.transpose(np.reshape(output['fx'],   (26,26)))
V = np.transpose(np.reshape(output['fz'],   (26,26)))
for i in range(26):
    for j in range(26):
        if U[i][j] < threshold:
            U[i][j] = 0
        if V[i][j] < threshold:
            V[i][j] = 0

speed = np.sqrt(U*U + V*V)


#speed[speed < 1] = 1
#speed = np.log10(speed)

lw = 4*speed/speed.max()


###########################
## Plotting ###############
###########################

# Make plot with vertical (default) colorbar
fig, ax = plt.subplots()


cax = plt.streamplot(X[:,0], Y[0,:], U, V, color=speed, linewidth=lw, cmap=cm.jet, arrowstyle="fancy")

xlabel = ax.set_xlabel('X (nm)', fontsize=labelSize)
ylabel = ax.set_ylabel('Z (nm)', fontsize=labelSize)
ax.tick_params(labelsize=axisSize)
ratio = (Y.max() - Y.min()) / (X.max() - X.min())
ax.set_aspect('equal')

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
#divider = make_axes_locatable(ax)
#caax = divider.append_axes("top", size="5%", pad=0.15)
cbaxes = fig.add_axes([0.125, 0.9, 0.775, 0.03])


#cbar = fig.colorbar(cax, cax=caax)
cbar = plt.colorbar(orientation='horizontal',cax=cbaxes)
#cbar.set_label('mV', fontsize=24)# vertically oriented colorbar
cbar_labels = np.arange(0,speed.max(),cbarInt)
cbar.ax.xaxis.set_ticks_position('top')
cbar.set_ticks(cbar_labels.astype(int))
cbar.ax.set_xticklabels(cbar_labels.astype(int),fontsize=axisSize)# vertically oriented colorbar

# set tick size
cbar.ax.tick_params('both', length=5, width=1.5, which='major')
# set axis width
for axis in ['top','bottom','left','right']:
  cbar.ax.spines[axis].set_linewidth(1.5)
#cbar.ax.invert_yaxis()

#cbar.ax.set_yticklabels([data_min, data_max], fontsize=24)# vertically oriented colorbar


plt.savefig(outputName, bbox_inches='tight',transparent=True)
plt.show()


