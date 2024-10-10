import os
import matplotlib.pyplot as plt
import numpy as np
from math import exp, log, sqrt
from scipy import interpolate
#--------------- data input -----------------------#
y0 = np.genfromtxt('Data/plot_1.tab',dtype=float, comments='#',unpack=True)
#--------------- writing plots -----------------#
#fig, ax = plt.subplots(figsize=(6,3.8))
#plt.subplots_adjust(bottom=0.20,right=0.75,left=0.15)
plt.title('e,E ->m,M')
plt.xlabel('$p_{cm}$',fontsize=16)
plt.ylabel('$v \sigma$[pb]',fontsize=16)
plt.xscale('linear')
plt.yscale('log')
xmin=1.000000E+01
xmax=2.000000E+02
plt.xlim(xmin,xmax)
plt.ylim(1.433872E+00,2.854165E+03)
xfirst=xmin+0.5*(xmax-xmin)/101
xlast=xmax -0.5*(xmax-xmin)/101
x0=np.linspace(xfirst,xlast,101)
plt.plot(x0, y0[0:101],alpha=0.7, label='$v \sigma$[pb]', linestyle='-')
plt.savefig('plot_1.pdf')
plt.show()

