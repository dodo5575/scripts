#!/usr/bin/env python3
# Usage python3 areaPerLipid.py pdb xst normal timestep lipidResname excludeVolume output
# author: Chen-Yu Li cli56@illinois.edu

import sys, optparse, json, os, glob, re, random, string, prody, scipy
import numpy as np
import pandas as pd
import seaborn as sb
import matplotlib.pyplot as plt
from subprocess import Popen, PIPE
from math import *
from scipy import stats
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import cm


def columns(x):
    return {
        'X': np.array([0, 5, 9]),
        'Y': np.array([0, 1, 9]),
        'Z': np.array([0, 1, 5]),
    }.get(x, 9)    # 9 is default if x not found


atoms = prody.parsePDB(sys.argv[1])

col = columns(sys.argv[3].upper())

xst = np.loadtxt(sys.argv[2], skiprows=2, usecols=col)

time = xst[:,0] * float(sys.argv[4]) * 10**(-6)

area = xst[:,1] * xst[:,2] - float(sys.argv[6])

#print(area)

lipidSel = 'resname %s and name P' % (sys.argv[5])

nLipid = atoms.select(lipidSel).numAtoms() / 2.0

#print(nLipid)

with open(sys.argv[7], 'a') as out:

    for i in range(len(time)):
        out.write('%f\t%f\n' % (time[i], area[i] / nLipid))



