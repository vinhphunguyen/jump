import numpy as np
from scipy import math
from numpy.random import random
from scipy.spatial import cKDTree
import matplotlib as mpl
import matplotlib.pyplot as plt
import pylab as pylab
import glob
import matplotlib.ticker as mtick
import matplotlib.colors as mcol
import time
import argparse
import sys, os

def set_size(width, fraction=1):
    """ Set aesthetic figure dimensions to avoid scaling in latex.

    Parameters
    ----------
    width: float
            Width in pts
    fraction: float
            Fraction of the width which you wish the figure to occupy

    Returns
    -------
    fig_dim: tuple
            Dimensions of figure in inches
    """

    if width == 'elsevier':
        width_pt = 468.
    elif width == 'beamer':
        width_pt = 307.28987
    else:
        width_pt = width

    # Width of figure
    fig_width_pt = width_pt * fraction

    # Convert from pt to inches
    inches_per_pt = 1 / 72.27

    # Golden ratio to set aesthetic figure height
    golden_ratio = 0.75 # (5**.5 - 1) / 2

    # Figure width in inches
    fig_width_in = fig_width_pt * inches_per_pt
    # Figure height in inches
    fig_height_in = fig_width_in * golden_ratio

    fig_dim = (fig_width_in, fig_height_in)

    return fig_dim

parser = argparse.ArgumentParser(description='Plot the evolution of variables.',
                                 add_help=False,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)


parser.add_argument('--quiet', '-q', action='store_const', const=True, default=False,
                    help='do not show plot.')
args = parser.parse_args()

# Using seaborn's style
plt.style.use('seaborn-colorblind') # dark_background,
width = 345

nice_fonts = {
        # Use LaTeX to write all text
        "text.usetex": True,
        "font.family": "serif",
        # Use 10pt font in plots, to match 10pt font in document
        "axes.labelsize": 10,
        "font.size": 10,
        # Make the legend/label fonts a little smaller
        "legend.fontsize": 8,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
}

mpl.rcParams.update(nice_fonts)

markersize = 3

#-------- Setup colours
#bgc = [(252./255.)**0.05, (141./255.)**0.05,(98./255.)**0.05]
#colors = brewer2mpl.get_map('Set2', 'qualitative', 8).mpl_colors


#-------------------------------------------------------------
# Dimer 5: Read Data

fnames = ['forces.txt']

dat = [pylab.genfromtxt(f,skip_header=1) for f in fnames]

x=[]
y=[]


fig, ax = plt.subplots(1, 1, figsize=set_size('elsevier', fraction=0.6))

lw = 1.5
ymax = 0.
xmax = 0.

for i, f in enumerate(fnames):
    x  = dat[i][:,0]
    fx = dat[i][:,1]*2
    fy = dat[i][:,2]*2
    x = x * 1e3 #- 6.07
    ax.plot(x, -fx, color='red', label='GPIC,shear',linewidth=lw)#, linestyle='dotted')
    ax.plot(x, -fy, color='black', label='GPIC,normal',linewidth=lw)#, linestyle='dotted')

    # xmax = max(xmax, np.max(x))
    # ymax = max(ymax, np.max(y))
    #if '460E' in f:
    #    xmax1 = max(xmax1, np.max(x[i]))
    #elif '700E' in f:
    #    xmax2 = max(xmax2, np.max(x[i]))
    #elif '900E' in f:
    #    xmax3 = max(xmax3, np.max(x[i]))


ax.set_xlabel(r'time in microsecond')
ax.set_ylabel(r'penetration [mm]')
# ax.set_ylim(0.0, 20.)
# ax.set_xlim(0, 50.)

ax.legend(loc=0)
plt.tight_layout()
plt.savefig('./high-velo-penetration.pdf', format='pdf')

#command = 'cp ' + 'water-break-plot.pdf' + ' /Users/vingu/Dropbox/my-writtings/papers-dropbox/karamelo/figures/'
#os.system( command )

if args.quiet == False:
    plt.show()
