#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from pylab import *
from scipy.signal import hilbert
import os, sys, pickle

# Plotting Parameters Setting
symbols = ['s','o','^','v','<','>','+','x','D','d'] # Symbol
lps = [k+'-' for k in ['o','^','v','<','>','s','+','x','D','d']] # Line + Symbol
colors= ['b','r','g','c','m','y','k','w'] # Color
ms = 5
ew = 2
lw = 1.0
rcParams['figure.figsize'] = 2*3.34646,5*1.67323
rcParams['ps.useafm'] = True
plt.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})

rcParams['pdf.fonttype'] = 42
matplotlib.rc('xtick.major', size=6)
matplotlib.rc('xtick.minor', size=3)
matplotlib.rc('ytick.major', size=6)
matplotlib.rc('ytick.minor', size=3)
matplotlib.rc('lines', markeredgewidth=0.5*2)
matplotlib.rc('font', size=8.5*2)

fig = plt.figure()
ax = fig.add_subplot(121)
ax.set_xlabel('Energy (eV)')
ax.set_ylabel('PDOS (Arb. Unit)')
vacuum = 0.11545125 * 27.2114 #vacuum
fermi = -0.077241 * 27.2114 # eV
shift = fermi-vacuum

DATA1 = np.loadtxt('all_pdos.txt')
plot(DATA1[:,1],DATA1[:,0]+shift,color=colors[0],linestyle='-')
fill_betweenx(DATA1[:,0]+shift, 0, DATA1[:,1], where=DATA1[:,0]<=0, color=colors[0])

ax.set_ylim([-8,0])
#plt.xticks([-6.0,-4.0,-2.0,0,2.0])                                                                                
plt.yticks([-8,-6,-4,-2,0])
minorticks_on()

ax = fig.add_subplot(122)

DATA1 = np.loadtxt('Co_pdos.txt')
plot(DATA1[:,1],DATA1[:,0]+shift,color=colors[0],linestyle='-')
fill_betweenx(DATA1[:,0]+shift, 0, DATA1[:,1], where=DATA1[:,0]<=0, color=colors[0])

DATA1 = np.loadtxt('N_pdos.txt')
plot(DATA1[:,1],DATA1[:,0]+shift,color=colors[0],linestyle='-')
fill_betweenx(DATA1[:,0]+shift, 0, DATA1[:,1], where=DATA1[:,0]<=0, color=colors[1])

DATA1 = np.loadtxt('O_pdos.txt')
plot(DATA1[:,1],DATA1[:,0]+shift,color=colors[0],linestyle='-')
fill_betweenx(DATA1[:,0]+shift, 0, DATA1[:,1], where=DATA1[:,0]<=0, color=colors[2])


ax.set_ylim([-8,0])
#plt.xticks([-6.0,-4.0,-2.0,0,2.0])
plt.yticks([-8,-6,-4,-2,0])
minorticks_on()

left  = 0.125  # the left side of the subplots of the figure
right = 0.95    # the right side of the subplots of the figure
bottom = 0.2   # the bottom of the subplots of the figure
top = 0.95      # the top of the subplots of the figure
wspace = 0.35   # the amount of width reserved for blank space between subplots
hspace = 0.35   # the amount of height reserved for white space between subplots
subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)

plt.savefig('Figure.pdf', format='pdf')
plt.show()
