#!/usr/bin/env python
# This script plots bulk ion concentration versus area.
# Usage: python plot_area_Mg_blockAverage.py 
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

titleSize = 20
labelSize = 16
axisSize = 12
xTickInt = 200
yTickInt = 0.05 


## Load Input File -------------------------------------------------
data0 = np.loadtxt('../Dropbox/Meeting_ppt/2014_6_2_meeting/square2plate_1MKCl_remove_Mg_npt40_Mg_area_bulkIon.dat')
data1 = np.loadtxt('/Users/chester/Dropbox/Meeting_ppt/2014_6_2_meeting/square2plate_1MKCl_reduce_Mg3_npt_extend25_mghh25_Mg_area_bulkIon.dat')
data2 = np.loadtxt('/Users/chester/Dropbox/Meeting_ppt/2014_6_2_meeting/square2plate_1MKCl_60_Mg_npt16_mghh25_Mg_area_bulkIon.dat')
data3 = np.loadtxt('/Users/chester/Dropbox/Meeting_ppt/2014_6_2_meeting/square2plate_1MKCl_73_Mg_npt15_mghh25_Mg_area_bulkIon.dat')
data4 = np.loadtxt('/Users/chester/Dropbox/Meeting_ppt/2014_6_2_meeting/square2plate_1MKCl_86_Mg_npt17_mghh25_Mg_area_bulkIon.dat')
data5 = np.loadtxt('/Users/chester/Dropbox/Meeting_ppt/2014_6_2_meeting/square2plate_1MKCl_99_Mg_npt17_mghh25_Mg_area_bulkIon.dat')
data6 = np.loadtxt('/Users/chester/Dropbox/Meeting_ppt/2014_6_2_meeting/square2plate_1MKCl_112_Mg_npt17_mghh25_Mg_area_bulkIon.dat')
data7 = np.loadtxt('/Users/chester/Dropbox/Meeting_ppt/2014_6_2_meeting/square2plate_1MKCl_npt_extend38_mghh25_Mg_area_bulkIon.dat')


## Set Output Name -------------------------------------------------

#blockSize = int(sys.argv[1])

outputName='Area_Mg_scatter.pdf'


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


## Pearson Correlation --------------------------------------------

#print scipy.stats.pearsonr(area_ba,P_ba)
#r0, p0 = scipy.stats.pearsonr(area0_ba,Mg0_ba)
r1, p1 = scipy.stats.pearsonr(area1_ba,Mg1_ba)
r2, p2 = scipy.stats.pearsonr(area2_ba,Mg2_ba)
r3, p3 = scipy.stats.pearsonr(area3_ba,Mg3_ba)
r4, p4 = scipy.stats.pearsonr(area4_ba,Mg4_ba)
r5, p5 = scipy.stats.pearsonr(area5_ba,Mg5_ba)
r6, p6 = scipy.stats.pearsonr(area6_ba,Mg6_ba)
r7, p7 = scipy.stats.pearsonr(area7_ba,Mg7_ba)


## Plotting -------------------------------------------------------

fig, ax = plt.subplots() 
plt.plot(area0_ba, Mg0_ba, 'v',color=(0.9961,0.5,0) ,label='SQ2 0 Mg$^{2+}$')
plt.plot(area1_ba, Mg1_ba, 'bo', label='SQ2 47 Mg$^{2+}$')
plt.plot(area2_ba, Mg2_ba, 'gv', label='SQ2 60 Mg$^{2+}$')
plt.plot(area3_ba, Mg3_ba, 'rh', label='SQ2 73 Mg$^{2+}$')
plt.plot(area4_ba, Mg4_ba, 'cd', label='SQ2 86 Mg$^{2+}$')
plt.plot(area5_ba, Mg5_ba, 'mp', label='SQ4 99 Mg$^{2+}$')
plt.plot(area6_ba, Mg6_ba, 'y*', label='SQ4 112 Mg$^{2+}$')
plt.plot(area7_ba, Mg7_ba, 'ks', label='SQ4 126 Mg$^{2+}$')
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

# Now add the legend with some customizations.
legend = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5),shadow=True,numpoints=1)
#legend = ax.legend(loc='upper right',shadow=True,numpoints=1)

# The frame is matplotlib.patches.Rectangle instance surrounding the legend.
#frame = legend.get_frame()
#frame.set_facecolor('0.90')

# Set the fontsize
for label in legend.get_texts():
    label.set_fontsize('large')

#for label in legend.get_lines():
#    label.set_linewidth(1.5)  # the legend line width



plt.savefig(outputName, bbox_inches='tight')
plt.show()




