#!/usr/bin/env python
# This script plots bulk ion concentration versus area.
# Usage: python plot_area_Mg_trace.py name xstFrequency timestep (blockSize) 
# Chen-Yu Li     cli56@illinois.edu
# 2014/5/29


import acor
import stats_blockAverage as sb
import math, random, sys, os, re, matplotlib 
import numpy as np
import matplotlib.pyplot as plt 
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


## get time ------------------------------------------------------

# interval (ns)
interval = float(sys.argv[2]) * float(sys.argv[3]) * 10**(-6)

time = []
for i in range(len(data3)):
    time.append(interval * i)


## block average -------------------------------------------------

if len(sys.argv) == 5:
    blockSize1 = int(sys.argv[4])
    blockSize2 = int(sys.argv[4])

    time_ba1 = sb.blockReduceFB(time,blockSize1)
    area_ba1 = sb.blockReduceFB(data3[:,0],blockSize1)
    time_ba2 = sb.blockReduceFB(time,blockSize2)
    Mg_ba2 = sb.blockReduceFB(data3[:,1],blockSize2)
    
elif len(sys.argv) == 4:

    acor_result1 = acor.acor(data3[:,0])
    blockSize1 = int(math.ceil(acor_result1[0]))
    acor_result2 = acor.acor(data3[:,1])
    blockSize2 = int(math.ceil(acor_result2[0]))
    
    time_ba1 = sb.blockReduceFB(time,blockSize1)
    area_ba1 = sb.blockReduceFB(data3[:,0],blockSize1)
    time_ba2 = sb.blockReduceFB(time,blockSize2)
    Mg_ba2 = sb.blockReduceFB(data3[:,1],blockSize2)


## Set Output Name -------------------------------------------------

s='_area_Mg_trace_b%db%d.pdf' % (blockSize1,blockSize2)
outputName = sys.argv[1]+s


## Plotting -------------------------------------------------------


plt.subplot(2, 1, 1)
plt.plot(time, data3[:,0], color=(0,0.9766,0.6016),alpha=0.5)
plt.plot(time_ba1,area_ba1,color=(0.0234,0.4101,0.8789),linewidth=2)
plt.title('Area', fontsize=titleSize)
plt.xlabel('Time (ns)', fontsize=labelSize)
plt.ylabel(r'Area ($\AA$)', fontsize=labelSize)
ax = plt.gca()
ax.tick_params(labelsize=axisSize)
ax.text(0.05, 0.95, 'ACT = %d' %(blockSize1),
        verticalalignment='top', horizontalalignment='left',
        transform=ax.transAxes,
        color='black', fontsize=axisSize)
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
ax.yaxis.grid()

plt.subplot(2, 1, 2)
plt.plot(time, data3[:,1],color=(0,0.9766,0.6016),alpha=0.5)
plt.plot(time_ba2,Mg_ba2,color=(0.0234,0.4101,0.8789),linewidth=2)
plt.title('Mg', fontsize=titleSize)
plt.xlabel('Time (ns)', fontsize=labelSize)
plt.ylabel('Conc. (M)', fontsize=labelSize)
ax = plt.gca()
ax.text(0.05, 0.95, 'ACT = %d' %(blockSize2),
        verticalalignment='top', horizontalalignment='left',
        transform=ax.transAxes,
        color='black', fontsize=axisSize)

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
ax.yaxis.grid()

#plt.subplot(2, 2, 2)
#plt.plot(time, data3[:,0], color=(0,0.9766,0.6016),alpha=0.5)
#plt.plot(time_ba2,area_ba2,color=(0.0234,0.4101,0.8789),linewidth=2)
#plt.title('Area', fontsize=titleSize)
#plt.xlabel('Time (ns)', fontsize=labelSize)
#plt.ylabel(r'Area ($\AA$)', fontsize=labelSize)
#ax = plt.gca()
#ax.tick_params(labelsize=axisSize)
#majorXLocator   = matplotlib.ticker.MultipleLocator(xTickInt*2)
#minorXLocator   = matplotlib.ticker.MultipleLocator(xTickInt)
#majorXFormatter = matplotlib.ticker.FormatStrFormatter('%d')
#ax.xaxis.set_major_locator(majorXLocator)
#ax.xaxis.set_minor_locator(minorXLocator)
#ax.xaxis.set_major_formatter(majorXFormatter)
##majorYLocator   = matplotlib.ticker.MultipleLocator(yTickInt*2)
##minorYLocator   = matplotlib.ticker.MultipleLocator(yTickInt)
##majorYFormatter = matplotlib.ticker.FormatStrFormatter('%d')
##ax.yaxis.set_major_locator(majorYLocator)
##ax.yaxis.set_minor_locator(minorYLocator)
##ax.yaxis.set_major_formatter(majorYFormatter)
#ax.yaxis.grid()
#
#plt.subplot(2, 2, 4)
#plt.plot(time, data3[:,1],color=(0,0.9766,0.6016),alpha=0.5)
#plt.plot(time_ba2,Mg_ba2,color=(0.0234,0.4101,0.8789),linewidth=2)
#plt.title('Mg', fontsize=titleSize)
#plt.xlabel('Time (ns)', fontsize=labelSize)
#plt.ylabel('Conc. (M)', fontsize=labelSize)
#ax = plt.gca()
#ax.tick_params(labelsize=axisSize)
#majorXLocator   = matplotlib.ticker.MultipleLocator(xTickInt*2)
#minorXLocator   = matplotlib.ticker.MultipleLocator(xTickInt)
#majorXFormatter = matplotlib.ticker.FormatStrFormatter('%d')
#ax.xaxis.set_major_locator(majorXLocator)
#ax.xaxis.set_minor_locator(minorXLocator)
#ax.xaxis.set_major_formatter(majorXFormatter)
##majorYLocator   = matplotlib.ticker.MultipleLocator(yTickInt*2)
##minorYLocator   = matplotlib.ticker.MultipleLocator(yTickInt)
##majorYFormatter = matplotlib.ticker.FormatStrFormatter('%d')
##ax.yaxis.set_major_locator(majorYLocator)
##ax.yaxis.set_minor_locator(minorYLocator)
##ax.yaxis.set_major_formatter(majorYFormatter)
#ax.yaxis.grid()


plt.savefig(outputName, bbox_inches='tight')
plt.show()




#!/usr/bin/env python
# This script plots bulk ion concentration versus area.
# Usage: python plot_area_Mg_trace.py name
