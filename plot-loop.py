import matplotlib.pyplot as plt
from matplotlib import rc,rcParams
rc('xtick', labelsize=14) 
rc('ytick', labelsize=14)
rc('font',**{'family':'serif','serif':['cmr10']})
# rcParams['mathtext.fontset'] = 'stix'
# rcParams['font.family'] = 'STIXGeneral'
rc('text', usetex=True)
# %matplotlib inline
from datetime import datetime
import sys
sys.path.append("/Users/andyreagan/work/2015/08-kitchentabletools/")
from dog.toys import *
import numpy as np

from matplotlib import patches

RdYlBu = np.genfromtxt("/Users/andyreagan/work/2014/11-foamLab-julia/RdYlBu.csv",delimiter=",")
def cmapper(cmap,val):
    return cmap[min(max(round(val*255),1),255)-1,:]

lookup_table = np.zeros([40000,4],dtype=np.int)
lookup_table[0:19999,:] = np.repeat([[0,1,3,2]],19999,axis=0)
# bottom right
lookup_table[19999:30000,:] = np.repeat([[0,1,2,3]],10001,axis=0)
# flip every 40th
lookup_table[19999:29999:40,:] = np.repeat([[0,1,3,2]],250,axis=0)
lookup_table[30000:40000,:] = np.repeat([[0,1,3,2]],10000,axis=0)
lookup_table[39,:] = [0,1,2,3]
lookup_table[39960:40000,:] = np.repeat([[0,1,2,3]],40,axis=0)
lookup_table[29960:29999,:] = np.repeat([[0,1,3,2]],39,axis=0)
def cell_point_order(cellID):
    # a lookup for how to reorder the point
    return lookup_table[cellID,:]

pts_x = np.genfromtxt("/Users/andyreagan/work/2014/11-foamLab-julia/cell_points_x.csv",delimiter=",")
pts_y = np.genfromtxt("/Users/andyreagan/work/2014/11-foamLab-julia/cell_points_y.csv",delimiter=",")

print(pts_x[0])
print(zip(pts_x[0],pts_y[0]))
a = np.array(zip(pts_x[0],pts_y[0]))
print(cell_point_order(0))
print(a[cell_point_order(0)])

def plot_loop_whole_save(T,fname,use_yticks=True,letter=""):
    # plot the normalized temp T
    # around the whole loop

    # 510.0 for a half-width figure
    textwidth = 510.0 # in pt
    textwidth_in = textwidth/72.27 # in in (on on)
    fig = plt.figure(figsize=(textwidth_in,textwidth_in))
    # width of the axes at .7 for a half-width figure works
    # need a 1/3 width figure, so take 2/3 of that.
    ax = fig.add_axes([0.2,0.2,0.7*2/3,0.7*2/3])
    
    for cellID in range(40000):
        pts = np.array(zip(pts_x[cellID],pts_y[cellID]))
        # println(pts)
        p = patches.Polygon(pts[cell_point_order(cellID)],closed=True,edgecolor="none",facecolor=cmapper(RdYlBu,T[cellID]),rasterized=False)
        ax.add_patch(p)
    ax.set_xlim(np.array([-1,1])*0.38)
    ax.set_ylim(np.array([-1,1])*0.38)

    if len(letter) > 0:
        props = dict(boxstyle='square', facecolor='white', alpha=1.0)
        ax.text(.34, -.34, letter,
                fontsize=16,
                verticalalignment='bottom',
                horizontalalignment='right',
                bbox=props)
    
    # ax.set_xlabel("Meters",fontsize=16)
    # if not use_yticks:
    #     ax.set_yticks([])
    # else:
    #     ax.set_ylabel("Meters",fontsize=16)

    ax.axis('off')
    
    # optionally save it to PDF:
    # mysavefig(fname)
    plt.savefig(fname,bbox_inches="tight",dpi=600)
    plt.close()

################################################################################
# do png. need rasterized=True and dpi=600.
    
# mode = np.genfromtxt("DMD-mode-2.csv")
# plot_loop_whole_save(mode,"DMD-mode-02-python-test-002.png")

# mode = np.genfromtxt("DMD-mode-21.csv")
# plot_loop_whole_save(mode,"DMD-mode-21-python-test-002.png",use_yticks=False)

# mode = np.genfromtxt("DMD-mode-79.csv")
# plot_loop_whole_save(mode,"DMD-mode-79-python-test-002.png")

################################################################################
# now for a test of SVG. using rasterized=False helped these actually save.

# mode = np.genfromtxt("DMD-mode-2.csv")
# plot_loop_whole_save(mode,"DMD-mode-02-python-test-2.svg")

# mode = np.genfromtxt("DMD-mode-21.csv")
# plot_loop_whole_save(mode,"DMD-mode-21-python-test-2.svg")

# mode = np.genfromtxt("DMD-mode-79.csv")
# plot_loop_whole_save(mode,"DMD-mode-79-python-test-2.svg")

################################################################################
# now for a test of PDF. using rasterized=False helped these actually save.

# mode = np.genfromtxt("DMD-mode-2.csv")
# plot_loop_whole_save(mode,"DMD-mode-02-python-test-003.pdf",letter="")

# mode = np.genfromtxt("DMD-mode-21.csv")
# plot_loop_whole_save(mode,"DMD-mode-21-python-test-003.pdf",letter="")

# mode = np.genfromtxt("DMD-mode-79.csv")
# plot_loop_whole_save(mode,"DMD-mode-79-python-test-003.pdf",letter="")

# mode = np.genfromtxt("DMD-mode-2.csv")
# plot_loop_whole_save(mode,"DMD-mode-02-python-test-004.pdf",letter="A")

# mode = np.genfromtxt("DMD-mode-21.csv")
# plot_loop_whole_save(mode,"DMD-mode-21-python-test-004.pdf",letter="B")

# mode = np.genfromtxt("DMD-mode-79.csv")
# plot_loop_whole_save(mode,"DMD-mode-79-python-test-004.pdf",letter="B")


mode = np.genfromtxt("DMD-mode-2.csv")
plot_loop_whole_save(mode,"DMD-mode-02-python-test-006.pdf",letter="A")

mode = np.genfromtxt("DMD-mode-21.csv")
plot_loop_whole_save(mode,"DMD-mode-21-python-test-006.pdf",letter="B")

mode = np.genfromtxt("DMD-mode-79.csv")
plot_loop_whole_save(mode,"DMD-mode-79-python-test-006.pdf",letter="B")
