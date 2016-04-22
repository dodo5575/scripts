#!/usr/bin/env python
# This script plots bulk ion concentration versus area.
# Usage: python plot_area_Mg_mean.py 
# Chen-Yu Li     cli56@illinois.edu
# 2014/5/30


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

titleSize = 28
labelSize = 24
axisSize = 20
xTickInt = 200
yTickInt = 0.05 
blockSize = 5000

## Load Input File -------------------------------------------------

data0 = np.loadtxt('../Dropbox/Meeting_ppt/2014_6_2_meeting/square2plate_1MKCl_remove_Mg_npt40_Mg_area_bulkIon.dat')
data1 = np.loadtxt('../Dropbox/Meeting_ppt/2014_6_2_meeting/square2plate_1MKCl_reduce_Mg3_npt_extend25_mghh25_Mg_area_bulkIon.dat')
data2 = np.loadtxt('../Dropbox/Meeting_ppt/2014_6_2_meeting/square2plate_1MKCl_60_Mg_npt16_mghh25_Mg_area_bulkIon.dat')
data3 = np.loadtxt('../Dropbox/Meeting_ppt/2014_6_2_meeting/square2plate_1MKCl_73_Mg_npt15_mghh25_Mg_area_bulkIon.dat')
data4 = np.loadtxt('../Dropbox/Meeting_ppt/2014_6_2_meeting/square2plate_1MKCl_86_Mg_npt17_mghh25_Mg_area_bulkIon.dat')
data5 = np.loadtxt('../Dropbox/Meeting_ppt/2014_6_2_meeting/square2plate_1MKCl_99_Mg_npt17_mghh25_Mg_area_bulkIon.dat')
data6 = np.loadtxt('../Dropbox/Meeting_ppt/2014_6_2_meeting/square2plate_1MKCl_112_Mg_npt17_mghh25_Mg_area_bulkIon.dat')
data7 = np.loadtxt('../Dropbox/Meeting_ppt/2014_6_2_meeting/square2plate_1MKCl_npt_extend38_mghh25_Mg_area_bulkIon.dat')


## block average -------------------------------------------------
acor_result0 = acor.acor(data0[:,0])
blockSize0 = int(math.ceil(acor_result0[0]))
area0_ba = sb.blockReduceFB(data0[:,0],blockSize0)
Mg0_ba = sb.blockReduceFB(data0[:,1],blockSize0)

acor_result1 = acor.acor(data1[:,0])
blockSize1 = int(math.ceil(acor_result1[0]))
area1_ba = sb.blockReduceFB(data1[:,0],blockSize1)
Mg1_ba = sb.blockReduceFB(data1[:,1],blockSize1)

acor_result2 = acor.acor(data2[:,0])
blockSize2 = int(math.ceil(acor_result2[0]))
area2_ba = sb.blockReduceFB(data2[:,0],blockSize2)
Mg2_ba = sb.blockReduceFB(data2[:,1],blockSize2)

acor_result3 = acor.acor(data3[:,0])
blockSize3 = int(math.ceil(acor_result3[0]))
area3_ba = sb.blockReduceFB(data3[:,0],blockSize3)
Mg3_ba = sb.blockReduceFB(data3[:,1],blockSize3)

acor_result4 = acor.acor(data4[:,0])
blockSize4 = int(math.ceil(acor_result4[0]))
area4_ba = sb.blockReduceFB(data4[:,0],blockSize4)
Mg4_ba = sb.blockReduceFB(data4[:,1],blockSize4)

acor_result5 = acor.acor(data5[:,0])
blockSize5 = int(math.ceil(acor_result5[0]))
area5_ba = sb.blockReduceFB(data5[:,0],blockSize5)
Mg5_ba = sb.blockReduceFB(data5[:,1],blockSize5)

acor_result6 = acor.acor(data6[:,0])
blockSize6 = int(math.ceil(acor_result6[0]))
area6_ba = sb.blockReduceFB(data6[:,0],blockSize6)
Mg6_ba = sb.blockReduceFB(data6[:,1],blockSize6)

acor_result7 = acor.acor(data7[:,0])
blockSize7 = int(math.ceil(acor_result7[0]))
area7_ba = sb.blockReduceFB(data7[:,0],blockSize7)
Mg7_ba = sb.blockReduceFB(data7[:,1],blockSize7)


## get mean and error ---------------------------------------------

x = [np.mean(data0[:,0]),np.mean(data1[:,0]),np.mean(data2[:,0]),np.mean(data3[:,0]),np.mean(data4[:,0]),np.mean(data5[:,0]),np.mean(data6[:,0]),np.mean(data7[:,0]),]

xErr = [scipy.stats.sem(area0_ba),scipy.stats.sem(area1_ba),scipy.stats.sem(area2_ba),scipy.stats.sem(area3_ba),scipy.stats.sem(area4_ba),scipy.stats.sem(area5_ba),scipy.stats.sem(area6_ba),scipy.stats.sem(area7_ba)]

y = [np.mean(data0[:,1]),np.mean(data1[:,1]),np.mean(data2[:,1]),np.mean(data3[:,1]),np.mean(data4[:,1]),np.mean(data5[:,1]),np.mean(data6[:,1]),np.mean(data7[:,1])]

yErr = [scipy.stats.sem(Mg0_ba),scipy.stats.sem(Mg1_ba),scipy.stats.sem(Mg2_ba),scipy.stats.sem(Mg3_ba),scipy.stats.sem(Mg4_ba),scipy.stats.sem(Mg5_ba),scipy.stats.sem(Mg6_ba),scipy.stats.sem(Mg7_ba)]


## Pearson Correlation --------------------------------------------

#print scipy.stats.pearsonr(area_ba,P_ba)
r1, p1 = scipy.stats.pearsonr(x,y)


## Set Output Name -------------------------------------------------

outputName='Area_Mg_mean.pdf'



## Plotting -------------------------------------------------------

fig, ax = plt.subplots() 
plt.plot(x, y, 'o',markersize=10, color=(0.0234,0.4101,0.8789))
plt.errorbar(x, y, xerr=xErr,yerr=yErr,fmt='bo',markersize=10, ecolor=(0.0234,0.4101,0.8789), capthick=3, capsize=6, elinewidth=2)

#plt.plot(area2_ba, Mg2_ba, 'go', label='SQ2 60 Mg$^{2+}$')
#plt.plot(area3_ba, Mg3_ba, 'ro', label='SQ2 73 Mg$^{2+}$')
#plt.plot(area4_ba, Mg4_ba, 'co', label='SQ2 86 Mg$^{2+}$')
#plt.plot(area5_ba, Mg5_ba, 'mo', label='SQ4 99 Mg$^{2+}$')
#plt.plot(area6_ba, Mg6_ba, 'yo', label='SQ4 112 Mg$^{2+}$')
#plt.plot(area7_ba, Mg7_ba, 'ko', label='SQ4 126 Mg$^{2+}$')
#plt.title('P', fontsize=titleSize)
plt.xlabel(r'Area ($\AA$)', fontsize=labelSize)
plt.ylabel('Conc. (M)', fontsize=labelSize)
#ax = plt.gca()
ax.tick_params(labelsize=axisSize)
majorXLocator   = matplotlib.ticker.MultipleLocator(xTickInt*2)
minorXLocator   = matplotlib.ticker.MultipleLocator(xTickInt)
majorXFormatter = matplotlib.ticker.FormatStrFormatter('%d')
ax.xaxis.set_major_locator(majorXLocator)
ax.xaxis.set_minor_locator(minorXLocator)
ax.xaxis.set_major_formatter(majorXFormatter)
majorYLocator   = matplotlib.ticker.MultipleLocator(yTickInt*2)
minorYLocator   = matplotlib.ticker.MultipleLocator(yTickInt)
majorYFormatter = matplotlib.ticker.FormatStrFormatter('%.2f')
ax.set_ylim([-0.05, 0.3])
ax.yaxis.set_major_locator(majorYLocator)
ax.yaxis.set_minor_locator(minorYLocator)
ax.yaxis.set_major_formatter(majorYFormatter)

ax.text(0.7, 0.9, 'R = %.2f' %(r1),
        verticalalignment='top', horizontalalignment='left',
        transform=ax.transAxes,
        color='black', fontsize=24)


# Now add the legend with some customizations.
#legend = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5),shadow=True,numpoints=1)

# The frame is matplotlib.patches.Rectangle instance surrounding the legend.
#frame = legend.get_frame()
#frame.set_facecolor('0.90')

# Set the fontsize
#for label in legend.get_texts():
#    label.set_fontsize('large')

#for label in legend.get_lines():
#    label.set_linewidth(1.5)  # the legend line width



plt.savefig(outputName, bbox_inches='tight')
plt.show()




