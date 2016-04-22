#!/usr/bin/env python
# This script plots bulk ion concentration versus area.
# Usage: python plot_area_bulkIon.py name blockSize
# Chen-Yu Li     cli56@illinois.edu
# 2014/5/29


import scipy
import acor
import stats_blockAverage as sb
import math, random, sys, os, re, matplotlib 
import numpy as np
import matplotlib.pyplot as plt 
from scipy import stats
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import cm
from numpy.random import randn


## Plotting Parameters ---------------------------------------------

titleSize = 20
labelSize = 16
axisSize = 12
xTickInt = 100
yTickInt = 100 


## Load Input File -------------------------------------------------

name1 = sys.argv[1]+'_P_area_bulkIon.dat'
name2 = sys.argv[1]+'_Cl_area_bulkIon.dat'
name3 = sys.argv[1]+'_Mg_area_bulkIon.dat'
name4 = sys.argv[1]+'_K_area_bulkIon.dat'


data1 = np.loadtxt(name1)
data2 = np.loadtxt(name2)
data3 = np.loadtxt(name3)
data4 = np.loadtxt(name4)


## Set Output Name -------------------------------------------------

blockSize = int(sys.argv[2])

s='_area_bulkIon_b%d.pdf' % (blockSize)
outputName = sys.argv[1]+s


## block average -------------------------------------------------

area_ba = sb.blockReduceFB(data1[:,0],blockSize)
P_ba = sb.blockReduceFB(data1[:,1],blockSize)
Cl_ba = sb.blockReduceFB(data2[:,1],blockSize)
Mg_ba = sb.blockReduceFB(data3[:,1],blockSize)
K_ba = sb.blockReduceFB(data4[:,1],blockSize)


## Pearson Correlation --------------------------------------------

#print scipy.stats.pearsonr(area_ba,P_ba)
Cl_r, Cl_p = scipy.stats.pearsonr(area_ba,Cl_ba)
Mg_r, Mg_p = scipy.stats.pearsonr(area_ba,Mg_ba)
K_r, K_p = scipy.stats.pearsonr(area_ba,K_ba)


## Plotting -------------------------------------------------------


plt.subplot(2, 2, 1)
plt.plot(area_ba, P_ba, '.')
plt.title('P', fontsize=titleSize)
plt.xlabel(r'Area ($\AA$)', fontsize=labelSize)
plt.ylabel('Conc. (M)', fontsize=labelSize)
ax = plt.gca()
ax.tick_params(labelsize=axisSize)
majorXLocator   = matplotlib.ticker.MultipleLocator(xTickInt*2)
minorXLocator   = matplotlib.ticker.MultipleLocator(xTickInt)
majorXFormatter = matplotlib.ticker.FormatStrFormatter('%d')
ax.xaxis.set_major_locator(majorXLocator)
ax.xaxis.set_minor_locator(minorXLocator)
ax.xaxis.set_major_formatter(majorXFormatter)
#majorYLocator   = matplotlib.ticker.MultipleLocator(yTickInt*2)
#minorYLocator   = matplotlib.ticker.MultipleLocator(yTickInt)
#majorYFormatter = matplotlib.ticker.FormatStrFormatter('%d')
#ax.yaxis.set_major_locator(majorYLocator)
#ax.yaxis.set_minor_locator(minorYLocator)
#ax.yaxis.set_major_formatter(majorYFormatter)


plt.subplot(2, 2, 2)
plt.plot(area_ba, Cl_ba, '.')
plt.title('Cl', fontsize=titleSize)
plt.xlabel(r'Area ($\AA$)', fontsize=labelSize)
plt.ylabel('Conc. (M)', fontsize=labelSize)
ax = plt.gca()
ax.tick_params(labelsize=axisSize)
majorXLocator   = matplotlib.ticker.MultipleLocator(xTickInt*2)
minorXLocator   = matplotlib.ticker.MultipleLocator(xTickInt)
majorXFormatter = matplotlib.ticker.FormatStrFormatter('%d')
ax.xaxis.set_major_locator(majorXLocator)
ax.xaxis.set_minor_locator(minorXLocator)
ax.xaxis.set_major_formatter(majorXFormatter)
#majorYLocator   = matplotlib.ticker.MultipleLocator(yTickInt*2)
#minorYLocator   = matplotlib.ticker.MultipleLocator(yTickInt)
#majorYFormatter = matplotlib.ticker.FormatStrFormatter('%d')
#ax.yaxis.set_major_locator(majorYLocator)
#ax.yaxis.set_minor_locator(minorYLocator)
#ax.yaxis.set_major_formatter(majorYFormatter)
ax.text(0.1, 0.9, 'R = %.2f' %(Cl_r),
        verticalalignment='top', horizontalalignment='left',
        transform=ax.transAxes,
        color='black', fontsize=15)


plt.subplot(2, 2, 3)
plt.plot(area_ba, Mg_ba, '.')
plt.title('Mg', fontsize=titleSize)
plt.xlabel(r'Area ($\AA$)', fontsize=labelSize)
plt.ylabel('Conc. (M)', fontsize=labelSize)
ax = plt.gca()
ax.tick_params(labelsize=axisSize)
majorXLocator   = matplotlib.ticker.MultipleLocator(xTickInt*2)
minorXLocator   = matplotlib.ticker.MultipleLocator(xTickInt)
majorXFormatter = matplotlib.ticker.FormatStrFormatter('%d')
ax.xaxis.set_major_locator(majorXLocator)
ax.xaxis.set_minor_locator(minorXLocator)
ax.xaxis.set_major_formatter(majorXFormatter)
#majorYLocator   = matplotlib.ticker.MultipleLocator(yTickInt*2)
#minorYLocator   = matplotlib.ticker.MultipleLocator(yTickInt)
#majorYFormatter = matplotlib.ticker.FormatStrFormatter('%d')
#ax.yaxis.set_major_locator(majorYLocator)
#ax.yaxis.set_minor_locator(minorYLocator)
#ax.yaxis.set_major_formatter(majorYFormatter)
ax.text(0.6, 0.9, 'R = %.2f' %(Mg_r),
        verticalalignment='top', horizontalalignment='left',
        transform=ax.transAxes,
        color='black', fontsize=15)


plt.subplot(2, 2, 4)
plt.plot(area_ba, K_ba, '.')
plt.title('K', fontsize=titleSize)
plt.xlabel(r'Area ($\AA$)', fontsize=labelSize)
plt.ylabel('Conc. (M)', fontsize=labelSize)
ax = plt.gca()
ax.tick_params(labelsize=axisSize)
majorXLocator   = matplotlib.ticker.MultipleLocator(xTickInt*2)
minorXLocator   = matplotlib.ticker.MultipleLocator(xTickInt)
majorXFormatter = matplotlib.ticker.FormatStrFormatter('%d')
ax.xaxis.set_major_locator(majorXLocator)
ax.xaxis.set_minor_locator(minorXLocator)
ax.xaxis.set_major_formatter(majorXFormatter)
#majorYLocator   = matplotlib.ticker.MultipleLocator(yTickInt*2)
#minorYLocator   = matplotlib.ticker.MultipleLocator(yTickInt)
#majorYFormatter = matplotlib.ticker.FormatStrFormatter('%d')
#ax.yaxis.set_major_locator(majorYLocator)
#ax.yaxis.set_minor_locator(minorYLocator)
#ax.yaxis.set_major_formatter(majorYFormatter)
ax.text(0.1, 0.9, 'R = %.2f' %(K_r),
        verticalalignment='top', horizontalalignment='left',
        transform=ax.transAxes,
        color='black', fontsize=15)



plt.savefig(outputName, bbox_inches='tight')
plt.show()




