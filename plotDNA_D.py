#!/usr/bin/env python
# This script generates colormap by averaging the dx file from a given axis.
# Usage: python plotStreamline.py input 


import math, random, sys, os, re, matplotlib, scipy
#matplotlib.use('Agg')
import scipy.ndimage
import numpy as np
import pandas as pd
import seaborn.apionly as sns
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

def pnploy(nVert, vertX, vertY, x, y):

    c = 0
    for i in range(nVert):
        if i == 0:
            j = nVert - 1
        else:
            j = i - 1
         
        if ((vertY[i] > y) != (vertY[j] > y)) and (x < (vertX[j] - vertX[i]) * (y - vertY[i]) / (vertY[j] - vertY[i]) + vertX[i]):
            c += 1 

    return c % 2


## Plotting Parameters ---------------------------------------------

titleSize = 32
labelSize = 28
axisSize = 24
xTickInt = 1
yTickInt = 1
cbarLow1 = 0
cbarInt1 = 20
cbarHigh1 = 90
D_seq = np.linspace(cbarLow1, cbarHigh1, 9)
cbarLow2 = 0
cbarInt2 = 0.1
cbarHigh2 = 0.5
Dmin = 0
Dmax = 90
Fmin = 0
Fmax = 0.5

lBound = -1.25
hBound =  1.25
threshold = 0
threshold1 = 60
dim = 'y'
xDim = 26
zDim = 26

current_convert = 1.6 * 10**(-19) * 10**(9) * 10**(9)

pdf = 'DNA_D.pdf'
png = 'DNA_D.png'
pdf_1 = 'DNA_D_1.pdf'
png_1 = 'DNA_D_1.png'

ionListMap = np.array([[0, 'POT', 1], 
                       [1, 'MG', 2], 
                       [2, 'CLA', -1]])

current_palette = sns.color_palette("bright")

####################
## Load Currents  ##
####################

inFileName = 'POTFlux.dat'

infile = open(inFileName, 'r')
#theLines = infile.readlines()

for i, line in enumerate(infile):
    if i == 2:
        gridDim = removeBlanksFromList(re.split('\s',line[15:len(line)]))
        #print gridDim
    if i == 3:
        gridNum = removeBlanksFromList(re.split('\s',line[13:len(line)]))
        #print gridNum
    elif i > 3:
        break

for i in range(len(gridNum)):
    gridDim[i] = float(gridDim[i])
    gridNum[i] = int(gridNum[i])

infile.close()

current_raw = np.zeros([3, gridNum[0]*gridNum[1]*gridNum[2], 7])

for ion in ionListMap:

    inFileName = '%sFlux.dat' % (str(ion[1]))
    tmp = np.loadtxt(inFileName, skiprows=6)
    #print(tmp)

    current_raw[int(ion[0])] = tmp.copy()

#print(current_raw)

current = np.zeros([gridNum[0]*gridNum[1]*gridNum[2], 7])
current[:,0:3] = current_raw[0,:,0:3].copy()

for ion in ionListMap:

    current[:,4] = current[:,4] + (current_raw[ion[0],:,4] * current_convert * float(ion[2]))
    current[:,5] = current[:,5] + (current_raw[ion[0],:,5] * current_convert * float(ion[2]))
    current[:,6] = current[:,6] + (current_raw[ion[0],:,6] * current_convert * float(ion[2]))

#print(current)

current_df = pd.DataFrame(current, columns=['x', 'y', 'z', 'dens', 'fx', 'fy', 'fz'])
#print(current_df)

sec_c = current_df[(current_df.z < hBound) & (current_df.z > lBound)].groupby(['x','y']).mean()
#print(sec_c)
#sys.exit()    
#X = current_df[(current_df.x < hBound) & (current_df.x > lBound)].groupby('y').mean().index
#Y = current_df[(current_df.x < hBound) & (current_df.x > lBound)].groupby('z').mean().index
#print(Y)

#FX = np.transpose(np.reshape(sec_c.fy,   (gridNum[1],gridNum[2])))
#FY = np.transpose(np.reshape(sec_c.fz,   (gridNum[1],gridNum[2])))
#speed = np.sqrt(FX*FX + FY*FY)
#lw = 5*speed/speed.max()

C = np.reshape(sec_c.fz, (gridNum[0],gridNum[1]))
C_T = np.transpose(C)
C_T_F = np.flipud(C_T)


############################
## Load Density ############
############################

inFileName = 'DNAFlux.dat'

infile = open(inFileName, 'r')

for i, line in enumerate(infile):
    if i == 2:
        gridDim = removeBlanksFromList(re.split('\s',line[15:len(line)]))
        #print gridDim
    if i == 3:
        gridNum = removeBlanksFromList(re.split('\s',line[13:len(line)]))
        #print gridNum
    elif i > 3:
        break

for i in range(len(gridNum)):
    gridDim[i] = float(gridDim[i])
    gridNum[i] = int(gridNum[i])

infile.close()

regexp = r"(-?[0-9]+[.][0-9]+)\s+(-?[0-9]+[.][0-9]+)\s+(-?[0-9]+[.][0-9]+)\s+(-?[0-9]+[.][0-9]+)\s+(-?[0-9]+[.][0-9]+)\s+(-?[0-9]+[.][0-9]+)\s+(-?[0-9]+[.][0-9]+)"

density = np.fromregex(inFileName, regexp,
                [('x', np.float), ('y', np.float), ('z', np.float), ('dens', np.float), ('fx', np.float), ('fy', np.float), ('fz', np.float) ])

density_df = pd.DataFrame(density)

sec_d = density_df[(density_df.z < hBound) & (density_df.z > lBound)].groupby(['x','y']).mean()

X = density_df.groupby('x').mean().index
Y = density_df.groupby('y').mean().index
Y_expand, X_expand = np.meshgrid(Y,X)
#print(X)
#print(Y)
#print(X_expand)
#print(Y_expand)
#sys.exit()

D = np.reshape(sec_d.dens, (gridNum[0],gridNum[1]))
D_T = np.transpose(D)
D_T_F = np.flipud(D_T)

#D_T = np.transpose(np.reshape(output['dens'], (xDim,zDim)))
# convert to M
D_M = 1.66 * D
#D_T_M = 1.66 * D_T

D_M[D_M < threshold] = 0


########################
## Find peak ###########
########################
D_T_thresh = np.copy(D_T)
D_T_thresh[D_T_thresh < threshold1] = 0

labeled_image, number_of_objects = scipy.ndimage.label(D_T_thresh)
#print(labeled_image,number_of_objects)
peak_slices = scipy.ndimage.find_objects(labeled_image)

centroids = []

for peak_slice in peak_slices:
    dy,dx  = peak_slice
    y,x = dy.start, dx.start
    #print(dy, dy.start, dx, dx.start)
    ind = np.where(D_T_thresh[peak_slice] == D_T_thresh[peak_slice].max())
    #print(ind)

    cy = ind[0]
    cx = ind[1]
    for i in range(len(cx)):

        #print y,x,cy[i],cx[i]
        centroids.append((y+cy[i],x+cx[i]))

centroids = np.array([centroids[1],
                      centroids[0],
                      centroids[3],
                      centroids[5],
                      centroids[4],
                      centroids[2]])
print(centroids)
#for y,x in centroids:
#    print(D_T[y,x])


########################
## Point in Polygon ####
########################

nVert = 6
vertY = centroids[:,0]
vertX = centroids[:,1]

#y = 20
#x = 24 
#print(pnploy(nVert, vertX, vertY, x, y))

inside = np.array([])
for y in range(len(Y)):
    for x in range(len(X)):

        if pnploy(nVert, vertX, vertY, x, y):
            inside = np.append(inside, [y,x])

inside = np.reshape(inside, (len(inside)/2, 2))
print("Number of grids inside:", len(inside))


########################
## Current in Polygon ##
########################

C_total = np.sum(C_T)
#print(C_total * 0.25)

C_in = 0
for y, x in inside:
    C_in += C_T[y, x]

#print(C_in/C_total, (C_total - C_in) / C_total)
ratio = [C_in / C_total, (C_total - C_in) / C_total]


###########################
## Plotting ###############
###########################

# Make plot with vertical (default) colorbar
#fig, ax = plt.subplots(figsize=(2.223, 12.189))
fig, ax = plt.subplots()

outDat = "DNA_D_PeakCoor.dat"
out = open(outDat, "w")


############################################# Density ######################################################

#print(D_M.min(), D_M.max())
print(D.min(), D.max())
#cax = ax.imshow(D_T_M, vmin = 0, vmax = 1.5, extent=[X.min(),X.max(),X.min(),X.max()], interpolation='nearest', cmap=cm.jet)
cax1 = ax.imshow(D_T_F, vmin = Dmin, vmax = Dmax, extent=[X_expand.min(),X_expand.max(),Y_expand.min(),Y_expand.max()], interpolation='nearest', cmap=cm.Blues)
#cax1 = plt.contourf(X_expand, Y_expand, D_M, D_seq, vmin = Dmin, vmax = Dmax, extent=[X_expand.min(),X_expand.max(),Y_expand.min(),Y_expand.max()], interpolation='nearest', cmap=cm.jet)
#, alpha=0.5

peakCoors = np.array([])
for y,x in centroids:
    y_coor = Y[y]
    x_coor = X[x]
    out.write("%10.3f%10.3f\n" % (x_coor,y_coor))
    peakCoors = np.append(peakCoors, [x_coor,y_coor])
    ax.plot(x_coor,y_coor,'kx',ms=10, markeredgewidth = 2)
    print(D_T[y, x])

out.close()

peakCoors = np.resize(peakCoors, (len(peakCoors)/2,2))
print(peakCoors)

ax.plot([peakCoors[0][0], peakCoors[1][0]], [peakCoors[0][1], peakCoors[1][1]], color = 'k', linewidth = 2)
ax.plot([peakCoors[1][0], peakCoors[2][0]], [peakCoors[1][1], peakCoors[2][1]], color = 'k', linewidth = 2)
ax.plot([peakCoors[2][0], peakCoors[3][0]], [peakCoors[2][1], peakCoors[3][1]], color = 'k', linewidth = 2)
ax.plot([peakCoors[3][0], peakCoors[4][0]], [peakCoors[3][1], peakCoors[4][1]], color = 'k', linewidth = 2)
ax.plot([peakCoors[4][0], peakCoors[5][0]], [peakCoors[4][1], peakCoors[5][1]], color = 'k', linewidth = 2)
ax.plot([peakCoors[5][0], peakCoors[0][0]], [peakCoors[5][1], peakCoors[0][1]], color = 'k', linewidth = 2)

#for y, x in inside:
#    y_coor = Y[y]
#    x_coor = X[x]
#    ax.plot(x_coor,y_coor,'wx',ms=10, markeredgewidth = 2)


xlabel = ax.set_xlabel('X (nm)', fontsize=labelSize)
ylabel = ax.set_ylabel('Y (nm)', fontsize=labelSize)
ax.tick_params(labelsize=axisSize)
ax.set_aspect('equal')

# set tick size
ax.tick_params('both', length=10, width=3, which='major')
ax.tick_params('both', length=5, width=2, which='minor')
# set axis width
for axis in ['top','bottom','left','right']:
  ax.spines[axis].set_linewidth(4)

#ax.set_xlim(0, 9)
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


############################################## Flux ######################################################
#
#print(speed.min(), speed.max())
#cax2 = plt.streamplot(X, Y, FX, FY, color=speed, linewidth=lw, cmap=cm.jet, arrowstyle="fancy", arrowsize=2.0, density=1, norm=matplotlib.colors.Normalize(vmin=Fmin, vmax=Fmax))
#
###plt.quiver(X, X, U_T, V_T, speed)
##xlabel = ax.set_xlabel('X (nm)', fontsize=labelSize)
##ylabel = ax.set_ylabel('Y (nm)', fontsize=labelSize)
##ax.tick_params(labelsize=axisSize)
##ax.set_aspect('equal')
##
### set tick size
##ax.tick_params('both', length=10, width=1.5, which='major')
##ax.tick_params('both', length=5, width=1, which='minor')
### set axis width
##for axis in ['top','bottom','left','right']:
##  ax.spines[axis].set_linewidth(1.5)
##
##majorXLocator   = matplotlib.ticker.MultipleLocator(xTickInt*2)
##minorXLocator   = matplotlib.ticker.MultipleLocator(xTickInt)
##majorXFormatter = matplotlib.ticker.FormatStrFormatter('%d')
##ax.xaxis.set_major_locator(majorXLocator)
##ax.xaxis.set_minor_locator(minorXLocator)
##ax.xaxis.set_major_formatter(majorXFormatter)
##
##majorYLocator   = matplotlib.ticker.MultipleLocator(yTickInt*2)
##minorYLocator   = matplotlib.ticker.MultipleLocator(yTickInt)
##majorYFormatter = matplotlib.ticker.FormatStrFormatter('%d')
##ax.yaxis.set_major_locator(majorYLocator)
##ax.yaxis.set_minor_locator(minorYLocator)
##ax.yaxis.set_major_formatter(majorYFormatter)


######################### Color bar 1 #####################################################################

# Add colorbar, make sure to specify tick locations to match desired ticklabels
#cbar = fig.colorbar(cax, ticks=[data_min, data_max])

# create an axes on the right side of ax. The width of cax will be 5%
# of ax and the padding between cax and ax will be fixed at 0.05 inch.
divider = make_axes_locatable(ax)
caax1 = divider.append_axes("right", size=0.18, pad=0.1)
#cbaxes = fig.add_axes([0.125, 0.9, 0.775, 0.03])


cbar1 = plt.colorbar(cax1, cax=caax1)
#cbar = plt.colorbar(orientation='horizontal',cax=cbaxes)
#cbar.set_label('mV', fontsize=24)# vertically oriented colorbar
cbar_labels1 = np.arange(cbarLow1, cbarHigh1, cbarInt1)
#for i in range(len(cbar_labels1)):
#    cbar_labels1[i] = round(cbar_labels1[i],2)
#cbar_labels_text1 = []
#for i in range(len(cbar_labels1)):
#    #if i == (len(cbar_labels) - 1):
#    if i == 0:
#        #cbar_labels_text.append(str(cbar_labels[i].astype(int))+' (nm$^{-3}$)')
#        cbar_labels_text1.append(str(cbar_labels1[i]))
#    else:
#        #cbar_labels_text.append(str(cbar_labels[i].astype(int)))
#        cbar_labels_text1.append(str(cbar_labels1[i]))

#cbar.ax.xaxis.set_ticks_position('top')
cbar1.set_ticks(cbar_labels1)
#cbar1.ax.set_yticklabels(cbar_labels1,fontsize=axisSize)# vertically oriented colorbar
cbar1.ax.tick_params(labelsize=axisSize)

ax.text(1.07, 1, '(nm$^{-3}$)',
       verticalalignment='top', horizontalalignment='left',
       transform=ax.transAxes,
       color='k', fontsize = axisSize)

# set tick size
cbar1.ax.tick_params('both', length=10, width=3, which='major')
# set axis width
for axis in ['top','bottom','left','right']:
  cbar1.ax.spines[axis].set_linewidth(4)
#cbar.ax.invert_yaxis()

#cbar.ax.set_yticklabels([data_min, data_max], fontsize=24)# vertically oriented colorbar

########################## Color bar 2 ####################################################################
#
## Add colorbar, make sure to specify tick locations to match desired ticklabels
##cbar = fig.colorbar(cax, ticks=[data_min, data_max])
#
## create an axes on the right side of ax. The width of cax will be 5%
## of ax and the padding between cax and ax will be fixed at 0.05 inch.
##divider = make_axes_locatable(ax)
##caax2 = divider.append_axes("right", size="5%", pad=1)
#cbaxes = fig.add_axes([0.85, 0.1, 0.023, 0.8]) #xmin, ymin, dx, dy
#
#
##cbar = fig.colorbar(cax, cax=caax)
##cbar2 = plt.colorbar(cax2, cax=cbaxes)
#cbar2 = plt.colorbar(cax=cbaxes)
##cbar2.set_clim(vmin=Fmin, vmax=Fmax)
#
##cbar.set_label('mV', fontsize=24)# vertically oriented colorbar
#cbar_labels2 = np.arange(cbarLow2, cbarHigh2, cbarInt2)
##cbar2.ax.xaxis.set_ticks_position('top')
#cbar2.set_ticks(cbar_labels2)
##cbar2.ax.set_yticklabels(cbar_labels2,fontsize=axisSize)# vertically oriented colorbar
#cbar2.ax.tick_params(labelsize=axisSize)
#
#ax.text(1.34, 1, '(nA/\nnm$^{2}$)',
#       verticalalignment='top', horizontalalignment='left',
#       transform=ax.transAxes,
#       color='k', fontsize = axisSize)
#
## set tick size
#cbar2.ax.tick_params('both', length=10, width=3, which='major')
## set axis width
#for axis in ['top','bottom','left','right']:
#  cbar2.ax.spines[axis].set_linewidth(4)




plt.savefig(pdf, bbox_inches='tight',transparent=True)
plt.savefig(png, bbox_inches='tight',transparent=True)
#plt.show()


###########################
## Plotting ###############
###########################

# Make plot with vertical (default) colorbar
#fig, ax = plt.subplots(figsize=(2.223, 12.189))
fig, ax = plt.subplots()

#outDat = "DNA_D_PeakCoor.dat"
#out = open(outDat, "w")


############################################# Density ######################################################

#print(D_M.min(), D_M.max())
#print(D.min(), D.max())
#cax = ax.imshow(D_T_M, vmin = 0, vmax = 1.5, extent=[X.min(),X.max(),X.min(),X.max()], interpolation='nearest', cmap=cm.jet)
cax1 = ax.imshow(D_T_F, vmin = Dmin, vmax = Dmax, extent=[X_expand.min(),X_expand.max(),Y_expand.min(),Y_expand.max()], interpolation='nearest', cmap=cm.Blues)
#cax1 = plt.contourf(X_expand, Y_expand, D_M, D_seq, vmin = Dmin, vmax = Dmax, extent=[X_expand.min(),X_expand.max(),Y_expand.min(),Y_expand.max()], interpolation='nearest', cmap=cm.jet)
#, alpha=0.5

peakCoors = np.array([])
for y,x in centroids:
    y_coor = Y[y]
    x_coor = X[x]
    #out.write("%10.3f%10.3f\n" % (x_coor,y_coor))
    peakCoors = np.append(peakCoors, [x_coor,y_coor])
    ax.plot(x_coor,y_coor,'kx',ms=10, markeredgewidth = 2)
    #print(D_T[y, x])

#out.close()

peakCoors = np.resize(peakCoors, (len(peakCoors)/2,2))
#print(peakCoors)

ax.plot([peakCoors[0][0], peakCoors[1][0]], [peakCoors[0][1], peakCoors[1][1]], color = 'k', linewidth = 2)
ax.plot([peakCoors[1][0], peakCoors[2][0]], [peakCoors[1][1], peakCoors[2][1]], color = 'k', linewidth = 2)
ax.plot([peakCoors[2][0], peakCoors[3][0]], [peakCoors[2][1], peakCoors[3][1]], color = 'k', linewidth = 2)
ax.plot([peakCoors[3][0], peakCoors[4][0]], [peakCoors[3][1], peakCoors[4][1]], color = 'k', linewidth = 2)
ax.plot([peakCoors[4][0], peakCoors[5][0]], [peakCoors[4][1], peakCoors[5][1]], color = 'k', linewidth = 2)
ax.plot([peakCoors[5][0], peakCoors[0][0]], [peakCoors[5][1], peakCoors[0][1]], color = 'k', linewidth = 2)

for y, x in inside:
    y_coor = Y[y]
    x_coor = X[x]
    ax.plot(x_coor,y_coor,'x',ms=10, markeredgewidth = 2, color=current_palette[1])


xlabel = ax.set_xlabel('X (nm)', fontsize=labelSize)
ylabel = ax.set_ylabel('Y (nm)', fontsize=labelSize)
ax.tick_params(labelsize=axisSize)
ax.set_aspect('equal')

# set tick size
ax.tick_params('both', length=10, width=3, which='major')
ax.tick_params('both', length=5, width=2, which='minor')
# set axis width
for axis in ['top','bottom','left','right']:
  ax.spines[axis].set_linewidth(4)

#ax.set_xlim(0, 9)
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


############################################## Flux ######################################################
#
#print(speed.min(), speed.max())
#cax2 = plt.streamplot(X, Y, FX, FY, color=speed, linewidth=lw, cmap=cm.jet, arrowstyle="fancy", arrowsize=2.0, density=1, norm=matplotlib.colors.Normalize(vmin=Fmin, vmax=Fmax))
#
###plt.quiver(X, X, U_T, V_T, speed)
##xlabel = ax.set_xlabel('X (nm)', fontsize=labelSize)
##ylabel = ax.set_ylabel('Y (nm)', fontsize=labelSize)
##ax.tick_params(labelsize=axisSize)
##ax.set_aspect('equal')
##
### set tick size
##ax.tick_params('both', length=10, width=1.5, which='major')
##ax.tick_params('both', length=5, width=1, which='minor')
### set axis width
##for axis in ['top','bottom','left','right']:
##  ax.spines[axis].set_linewidth(1.5)
##
##majorXLocator   = matplotlib.ticker.MultipleLocator(xTickInt*2)
##minorXLocator   = matplotlib.ticker.MultipleLocator(xTickInt)
##majorXFormatter = matplotlib.ticker.FormatStrFormatter('%d')
##ax.xaxis.set_major_locator(majorXLocator)
##ax.xaxis.set_minor_locator(minorXLocator)
##ax.xaxis.set_major_formatter(majorXFormatter)
##
##majorYLocator   = matplotlib.ticker.MultipleLocator(yTickInt*2)
##minorYLocator   = matplotlib.ticker.MultipleLocator(yTickInt)
##majorYFormatter = matplotlib.ticker.FormatStrFormatter('%d')
##ax.yaxis.set_major_locator(majorYLocator)
##ax.yaxis.set_minor_locator(minorYLocator)
##ax.yaxis.set_major_formatter(majorYFormatter)


######################### Color bar 1 #####################################################################

# Add colorbar, make sure to specify tick locations to match desired ticklabels
#cbar = fig.colorbar(cax, ticks=[data_min, data_max])

# create an axes on the right side of ax. The width of cax will be 5%
# of ax and the padding between cax and ax will be fixed at 0.05 inch.
divider = make_axes_locatable(ax)
caax1 = divider.append_axes("right", size=0.18, pad=0.1)
#cbaxes = fig.add_axes([0.125, 0.9, 0.775, 0.03])


cbar1 = plt.colorbar(cax1, cax=caax1)
#cbar = plt.colorbar(orientation='horizontal',cax=cbaxes)
#cbar.set_label('mV', fontsize=24)# vertically oriented colorbar
cbar_labels1 = np.arange(cbarLow1, cbarHigh1, cbarInt1)
#for i in range(len(cbar_labels1)):
#    cbar_labels1[i] = round(cbar_labels1[i],2)
#cbar_labels_text1 = []
#for i in range(len(cbar_labels1)):
#    #if i == (len(cbar_labels) - 1):
#    if i == 0:
#        #cbar_labels_text.append(str(cbar_labels[i].astype(int))+' (nm$^{-3}$)')
#        cbar_labels_text1.append(str(cbar_labels1[i]))
#    else:
#        #cbar_labels_text.append(str(cbar_labels[i].astype(int)))
#        cbar_labels_text1.append(str(cbar_labels1[i]))

#cbar.ax.xaxis.set_ticks_position('top')
cbar1.set_ticks(cbar_labels1)
#cbar1.ax.set_yticklabels(cbar_labels1,fontsize=axisSize)# vertically oriented colorbar
cbar1.ax.tick_params(labelsize=axisSize)

ax.text(1.07, 1, '(nm$^{-3}$)',
       verticalalignment='top', horizontalalignment='left',
       transform=ax.transAxes,
       color='k', fontsize = axisSize)

# set tick size
cbar1.ax.tick_params('both', length=10, width=3, which='major')
# set axis width
for axis in ['top','bottom','left','right']:
  cbar1.ax.spines[axis].set_linewidth(4)
#cbar.ax.invert_yaxis()

#cbar.ax.set_yticklabels([data_min, data_max], fontsize=24)# vertically oriented colorbar

########################## Color bar 2 ####################################################################
#
## Add colorbar, make sure to specify tick locations to match desired ticklabels
##cbar = fig.colorbar(cax, ticks=[data_min, data_max])
#
## create an axes on the right side of ax. The width of cax will be 5%
## of ax and the padding between cax and ax will be fixed at 0.05 inch.
##divider = make_axes_locatable(ax)
##caax2 = divider.append_axes("right", size="5%", pad=1)
#cbaxes = fig.add_axes([0.85, 0.1, 0.023, 0.8]) #xmin, ymin, dx, dy
#
#
##cbar = fig.colorbar(cax, cax=caax)
##cbar2 = plt.colorbar(cax2, cax=cbaxes)
#cbar2 = plt.colorbar(cax=cbaxes)
##cbar2.set_clim(vmin=Fmin, vmax=Fmax)
#
##cbar.set_label('mV', fontsize=24)# vertically oriented colorbar
#cbar_labels2 = np.arange(cbarLow2, cbarHigh2, cbarInt2)
##cbar2.ax.xaxis.set_ticks_position('top')
#cbar2.set_ticks(cbar_labels2)
##cbar2.ax.set_yticklabels(cbar_labels2,fontsize=axisSize)# vertically oriented colorbar
#cbar2.ax.tick_params(labelsize=axisSize)
#
#ax.text(1.34, 1, '(nA/\nnm$^{2}$)',
#       verticalalignment='top', horizontalalignment='left',
#       transform=ax.transAxes,
#       color='k', fontsize = axisSize)
#
## set tick size
#cbar2.ax.tick_params('both', length=10, width=3, which='major')
## set axis width
#for axis in ['top','bottom','left','right']:
#  cbar2.ax.spines[axis].set_linewidth(4)




plt.savefig(pdf_1, bbox_inches='tight',transparent=True)
plt.savefig(png_1, bbox_inches='tight',transparent=True)
#plt.show()


###########################
## Plotting ###############
###########################

## Plotting Parameters ---------------------------------------------

cbarLow1 = 0
cbarInt1 = 0.05
cbarHigh1 = 0.2
D_seq = np.linspace(cbarLow1, cbarHigh1, 5)
Dmin = 0
Dmax = 0.2

pdf_2 = 'DNA_D_2.pdf'
png_2 = 'DNA_D_2.png'


# Make plot with vertical (default) colorbar
#fig, ax = plt.subplots(figsize=(2.223, 12.189))
fig, ax = plt.subplots()

#outDat = "DNA_D_PeakCoor.dat"
#out = open(outDat, "w")


############################################# Density ######################################################

print(C_T_F.min(), C_T_F.max())
#print(D.min(), D.max())
#cax = ax.imshow(D_T_M, vmin = 0, vmax = 1.5, extent=[X.min(),X.max(),X.min(),X.max()], interpolation='nearest', cmap=cm.jet)
cax1 = ax.imshow(C_T_F, vmin = Dmin, vmax = Dmax, extent=[X_expand.min(),X_expand.max(),Y_expand.min(),Y_expand.max()], interpolation='nearest', cmap=cm.Reds)
#cax1 = plt.contourf(X_expand, Y_expand, D_M, D_seq, vmin = Dmin, vmax = Dmax, extent=[X_expand.min(),X_expand.max(),Y_expand.min(),Y_expand.max()], interpolation='nearest', cmap=cm.jet)
#, alpha=0.5

peakCoors = np.array([])
for y,x in centroids:
    y_coor = Y[y]
    x_coor = X[x]
    #out.write("%10.3f%10.3f\n" % (x_coor,y_coor))
    peakCoors = np.append(peakCoors, [x_coor,y_coor])
    ax.plot(x_coor,y_coor,'kx',ms=10, markeredgewidth = 2)
    #print(D_T[y, x])

#out.close()

peakCoors = np.resize(peakCoors, (len(peakCoors)/2,2))
#print(peakCoors)

ax.plot([peakCoors[0][0], peakCoors[1][0]], [peakCoors[0][1], peakCoors[1][1]], color = 'k', linewidth = 2)
ax.plot([peakCoors[1][0], peakCoors[2][0]], [peakCoors[1][1], peakCoors[2][1]], color = 'k', linewidth = 2)
ax.plot([peakCoors[2][0], peakCoors[3][0]], [peakCoors[2][1], peakCoors[3][1]], color = 'k', linewidth = 2)
ax.plot([peakCoors[3][0], peakCoors[4][0]], [peakCoors[3][1], peakCoors[4][1]], color = 'k', linewidth = 2)
ax.plot([peakCoors[4][0], peakCoors[5][0]], [peakCoors[4][1], peakCoors[5][1]], color = 'k', linewidth = 2)
ax.plot([peakCoors[5][0], peakCoors[0][0]], [peakCoors[5][1], peakCoors[0][1]], color = 'k', linewidth = 2)

#for y, x in inside:
#    y_coor = Y[y]
#    x_coor = X[x]
#    ax.plot(x_coor,y_coor,'gx',ms=10, markeredgewidth = 2)


xlabel = ax.set_xlabel('X (nm)', fontsize=labelSize)
ylabel = ax.set_ylabel('Y (nm)', fontsize=labelSize)
ax.tick_params(labelsize=axisSize)
ax.set_aspect('equal')

# set tick size
ax.tick_params('both', length=10, width=3, which='major')
ax.tick_params('both', length=5, width=2, which='minor')
# set axis width
for axis in ['top','bottom','left','right']:
  ax.spines[axis].set_linewidth(4)

#ax.set_xlim(0, 9)
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


############################################## Flux ######################################################
#
#print(speed.min(), speed.max())
#cax2 = plt.streamplot(X, Y, FX, FY, color=speed, linewidth=lw, cmap=cm.jet, arrowstyle="fancy", arrowsize=2.0, density=1, norm=matplotlib.colors.Normalize(vmin=Fmin, vmax=Fmax))
#
###plt.quiver(X, X, U_T, V_T, speed)
##xlabel = ax.set_xlabel('X (nm)', fontsize=labelSize)
##ylabel = ax.set_ylabel('Y (nm)', fontsize=labelSize)
##ax.tick_params(labelsize=axisSize)
##ax.set_aspect('equal')
##
### set tick size
##ax.tick_params('both', length=10, width=1.5, which='major')
##ax.tick_params('both', length=5, width=1, which='minor')
### set axis width
##for axis in ['top','bottom','left','right']:
##  ax.spines[axis].set_linewidth(1.5)
##
##majorXLocator   = matplotlib.ticker.MultipleLocator(xTickInt*2)
##minorXLocator   = matplotlib.ticker.MultipleLocator(xTickInt)
##majorXFormatter = matplotlib.ticker.FormatStrFormatter('%d')
##ax.xaxis.set_major_locator(majorXLocator)
##ax.xaxis.set_minor_locator(minorXLocator)
##ax.xaxis.set_major_formatter(majorXFormatter)
##
##majorYLocator   = matplotlib.ticker.MultipleLocator(yTickInt*2)
##minorYLocator   = matplotlib.ticker.MultipleLocator(yTickInt)
##majorYFormatter = matplotlib.ticker.FormatStrFormatter('%d')
##ax.yaxis.set_major_locator(majorYLocator)
##ax.yaxis.set_minor_locator(minorYLocator)
##ax.yaxis.set_major_formatter(majorYFormatter)


######################### Color bar 1 #####################################################################

# Add colorbar, make sure to specify tick locations to match desired ticklabels
#cbar = fig.colorbar(cax, ticks=[data_min, data_max])

# create an axes on the right side of ax. The width of cax will be 5%
# of ax and the padding between cax and ax will be fixed at 0.05 inch.
divider = make_axes_locatable(ax)
caax1 = divider.append_axes("right", size=0.18, pad=0.1)
#cbaxes = fig.add_axes([0.125, 0.9, 0.775, 0.03])


cbar1 = plt.colorbar(cax1, cax=caax1)
#cbar = plt.colorbar(orientation='horizontal',cax=cbaxes)
#cbar.set_label('mV', fontsize=24)# vertically oriented colorbar
cbar_labels1 = np.arange(cbarLow1, cbarHigh1, cbarInt1)
#for i in range(len(cbar_labels1)):
#    cbar_labels1[i] = round(cbar_labels1[i],2)
#cbar_labels_text1 = []
#for i in range(len(cbar_labels1)):
#    #if i == (len(cbar_labels) - 1):
#    if i == 0:
#        #cbar_labels_text.append(str(cbar_labels[i].astype(int))+' (nm$^{-3}$)')
#        cbar_labels_text1.append(str(cbar_labels1[i]))
#    else:
#        #cbar_labels_text.append(str(cbar_labels[i].astype(int)))
#        cbar_labels_text1.append(str(cbar_labels1[i]))

#cbar.ax.xaxis.set_ticks_position('top')
cbar1.set_ticks(cbar_labels1)
#cbar1.ax.set_yticklabels(cbar_labels1,fontsize=axisSize)# vertically oriented colorbar
cbar1.ax.tick_params(labelsize=axisSize)

ax.text(1.07, 1, '(nA /\nnm$^{2}$)',
       verticalalignment='top', horizontalalignment='left',
       transform=ax.transAxes,
       color='k', fontsize = axisSize)

# set tick size
cbar1.ax.tick_params('both', length=10, width=3, which='major')
# set axis width
for axis in ['top','bottom','left','right']:
  cbar1.ax.spines[axis].set_linewidth(4)
#cbar.ax.invert_yaxis()

#cbar.ax.set_yticklabels([data_min, data_max], fontsize=24)# vertically oriented colorbar

########################## Color bar 2 ####################################################################
#
## Add colorbar, make sure to specify tick locations to match desired ticklabels
##cbar = fig.colorbar(cax, ticks=[data_min, data_max])
#
## create an axes on the right side of ax. The width of cax will be 5%
## of ax and the padding between cax and ax will be fixed at 0.05 inch.
##divider = make_axes_locatable(ax)
##caax2 = divider.append_axes("right", size="5%", pad=1)
#cbaxes = fig.add_axes([0.85, 0.1, 0.023, 0.8]) #xmin, ymin, dx, dy
#
#
##cbar = fig.colorbar(cax, cax=caax)
##cbar2 = plt.colorbar(cax2, cax=cbaxes)
#cbar2 = plt.colorbar(cax=cbaxes)
##cbar2.set_clim(vmin=Fmin, vmax=Fmax)
#
##cbar.set_label('mV', fontsize=24)# vertically oriented colorbar
#cbar_labels2 = np.arange(cbarLow2, cbarHigh2, cbarInt2)
##cbar2.ax.xaxis.set_ticks_position('top')
#cbar2.set_ticks(cbar_labels2)
##cbar2.ax.set_yticklabels(cbar_labels2,fontsize=axisSize)# vertically oriented colorbar
#cbar2.ax.tick_params(labelsize=axisSize)
#
#ax.text(1.34, 1, '(nA/\nnm$^{2}$)',
#       verticalalignment='top', horizontalalignment='left',
#       transform=ax.transAxes,
#       color='k', fontsize = axisSize)
#
## set tick size
#cbar2.ax.tick_params('both', length=10, width=3, which='major')
## set axis width
#for axis in ['top','bottom','left','right']:
#  cbar2.ax.spines[axis].set_linewidth(4)




plt.savefig(pdf_2, bbox_inches='tight',transparent=True)
plt.savefig(png_2, bbox_inches='tight',transparent=True)
#plt.show()


###########################
## Plotting ###############
###########################

out = open("ratio.dat", "w")
out.write("%f\t%f\n" % (ratio[0], ratio[1]))
out.close

## Plotting Parameters ---------------------------------------------

current_palette = sns.color_palette("pastel")

# The slices will be ordered and plotted counter-clockwise.
labels = ['in', 'out']
explode = [0, 0] # only "explode" the 2nd slice (i.e. 'Hogs')

pdf_3 = 'DNA_D_3.pdf'
png_3 = 'DNA_D_3.png'


fig, ax = plt.subplots()

patches, texts, autotexts = plt.pie(ratio, explode=explode, labels=labels, colors=current_palette[0:2],
textprops={'fontsize': labelSize},
wedgeprops={'linewidth' : 3 },
                autopct='%1.1f%%', shadow=True, startangle=90)

#print(texts)
for i in range(len(texts)):
    texts[i].set_fontsize(labelSize)

ax.set_aspect('equal')

plt.savefig(pdf_3, bbox_inches='tight',transparent=True)
plt.savefig(png_3, bbox_inches='tight',transparent=True)
plt.show()


