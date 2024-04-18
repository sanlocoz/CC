# -*- coding: utf-8 -*-
"""
Created on Mon Apr 15 09:53:24 2024

@author: user
"""
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 08:49:24 2024

@author: user
"""
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  4 17:46:34 2024

@author: user
"""
# %% Import library 
import numpy as np
import seaborn as sns
import numpy as np
import os
import re
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import scipy.stats as stats
import warnings
import pdb


# %% Import the dataset of measurements
dir = "D:\Wageningen\Period 5\Catchment\Module 3\Practical\WALRUS\data\\viz"
os.chdir(dir)
files = os.listdir(dir)
output_dir = "D:\Wageningen\Period 5\Catchment\Module 3\Practical\WALRUS\output\fig"

forcing = pd.read_csv(files[3], index_col= 0, delim_whitespace=True).iloc[365*24:,]

forcing['DateTime'] = forcing['date'].apply(lambda x: pd.to_datetime(str(x), format='%Y%m%d%H'))

QR = pd.read_csv(files[4], index_col= 0, delim_whitespace=True).iloc[365*24:,]
QWC = pd.read_csv(files[5], index_col= 0, delim_whitespace=True).iloc[365*24:-1,]
QWE = pd.read_csv(files[6], index_col= 0, delim_whitespace=True).iloc[365*24:-1,]
QWQ = pd.read_csv(files[7], index_col= 0, delim_whitespace=True).iloc[365*24:-1,]

dR = pd.read_csv(files[8], index_col= 0, delim_whitespace=True).iloc[365*24:,]
dC = pd.read_csv(files[0], index_col= 0, delim_whitespace=True).iloc[365*24:-1,]
dE = pd.read_csv(files[1], index_col= 0, delim_whitespace=True).iloc[365*24:-1,]
dQ = pd.read_csv(files[2], index_col= 0, delim_whitespace=True).iloc[365*24:-1,]

#%%
def my_plotter(ax, data1, data2, param_dict):
    """
    A helper function to make a graph.
    """
    out = ax.plot(data1, data2, **param_dict)
    return out

#%%
## PLOT FORCING
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

fig, ax = plt.subplots(2, 1, sharex= True, figsize = (14,12))

my_plotter(ax[0], forcing["DateTime"], forcing["P"], {"color" : "blue", "label":"Precipitation"})
cdf = mpl.dates.ConciseDateFormatter(ax[0].xaxis.get_major_locator())
ax[0].xaxis.set_major_formatter(cdf)
ax[0].set_ylim([17.5*5, 0])

ax3 = ax[0].twinx()

inset_ax = ax3.inset_axes(
   [0.1, 0.3, 0.15, 0.6],  # [x, y, width, height] w.r.t. axes
    xlim=[forcing["DateTime"].iloc[4000], forcing["DateTime"].iloc[4500]], ylim=[0, 0.15], # sets viewport & tells relation to main axes
    xticklabels=[], yticklabels=[]
)

inset2_ax = ax3.inset_axes(
   [0.75, 0.25, 0.2, 0.4],  # [x, y, width, height] w.r.t. axes
    xlim=[forcing["DateTime"].iloc[7000], forcing["DateTime"].iloc[8000]], ylim=[0.06, 0.2], # sets viewport & tells relation to main axes
    xticklabels=[], yticklabels=[]
)

ax3.indicate_inset_zoom(inset_ax, edgecolor="blue")
ax3.indicate_inset_zoom(inset2_ax, edgecolor="blue")

for axx in ax3, inset_ax, inset2_ax:
    my_plotter(axx, forcing["DateTime"], forcing["Q"], {"color" : "black", "label":"Q observation", "linewidth":1})
    ax3.set_ylim([0,0.59784*1.7])

    my_plotter(axx, forcing["DateTime"], QR["d50"], {"color" : "red", "label":"Q recession, KGE = 0.75", "alpha":0.6})
    plt.fill_between(forcing["DateTime"], QR["d25"], QR["d75"], alpha=0.3, color = "red")

    my_plotter(axx, forcing["DateTime"], QWQ["d50"], {"color" : "green", "label":"Q WALRUS, Q-driven, KGE = 0.68", "alpha":0.6})
    plt.fill_between(forcing["DateTime"], QWQ["d25"], QWQ["d75"], alpha=0.3, color = "green")
    
    my_plotter(axx, forcing["DateTime"], QWE["d50"], {"color" : "orange", "label":"Q WALRUS, ET-driven, KGE = 0.44", "alpha":0.6})
    plt.fill_between(forcing["DateTime"], QWE["d25"], QWE["d75"], alpha=0.3, color = "orange")
    
    my_plotter(axx, forcing["DateTime"], QWC["d50"], {"color" : "blue", "label":"Q WALRUS, Composite, KGE = 0.70", "alpha":0.6})
    plt.fill_between(forcing["DateTime"], QWC["d25"], QWC["d75"], alpha=0.3, color = "blue")

ax[0].legend(loc = 2)
ax3.legend(loc = 1)


from math import isclose; import matplotlib.pyplot as plt
if not isclose(inset_ax._get_aspect_ratio(), ax3._get_aspect_ratio()):
    print("chosen inset x/ylim & width/height skew aspect w.r.t. main axes!")

ax[0].set_ylabel("Precipitation (mm/h)")

ax3.set_ylabel("Streamflow (mm/h)")

#%%


## PLOT FORCING
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

my_plotter(ax[1], forcing["DateTime"], forcing["ETpot"], {"color" : "blue", "label":"ET potential"})
cdf = mpl.dates.ConciseDateFormatter(ax[0].xaxis.get_major_locator())
ax[1].xaxis.set_major_formatter(cdf)
ax[1].set_ylim([5, 0])

ax31 = ax[1].twinx()

# inset_ax = ax31.inset_axes(
#    [0.1, 0.3, 0.15, 0.6],  # [x, y, width, height] w.r.t. axes
#     xlim=[forcing["DateTime"].iloc[4000], forcing["DateTime"].iloc[4500]], ylim=[0, 0.15], # sets viewport & tells relation to main axes
#     xticklabels=[], yticklabels=[]
# )

# inset2_ax = ax31.inset_axes(
#    [0.75, 0.25, 0.2, 0.4],  # [x, y, width, height] w.r.t. axes
#     xlim=[forcing["DateTime"].iloc[7000], forcing["DateTime"].iloc[8000]], ylim=[0.06, 0.2], # sets viewport & tells relation to main axes
#     xticklabels=[], yticklabels=[]
# )

# ax31.indicate_inset_zoom(inset_ax, edgecolor="blue")
# ax31.indicate_inset_zoom(inset2_ax, edgecolor="blue")
constant = 150
for axx in ax31,:

    my_plotter(axx, forcing["DateTime"], dR["d50"]+constant, {"color" : "red", "label":"Storage recession", "alpha":1})
    plt.fill_between(forcing["DateTime"], dR["d25"]+constant, dR["d75"]+constant, alpha=0.2, color = "red")

    my_plotter(axx, forcing["DateTime"], dQ["d50"], {"color" : "green", "label":"Storage deficit WALRUS, Q-driven", "alpha":0.6})
    plt.fill_between(forcing["DateTime"], dQ["d25"], dQ["d75"], alpha=0.3, color = "green")
    
    my_plotter(axx, forcing["DateTime"], dE["d50"], {"color" : "orange", "label":"Storage deficit WALRUS, ET-driven", "alpha":0.6})
    plt.fill_between(forcing["DateTime"], dE["d25"], dE["d75"], alpha=0.3, color = "orange")
    
    my_plotter(axx, forcing["DateTime"], dC["d50"], {"color" : "blue", "label":"Storage deficit WALRUS, Composite", "alpha":0.6})
    plt.fill_between(forcing["DateTime"], dC["d25"], dC["d75"], alpha=0.3, color = "blue")

ax[1].legend(loc = 2)
ax31.legend(loc = 1)
ax31.set_ylim([0,700])

ax[1].set_ylabel("Evapotranspiration (mm/h)")

ax31.set_ylabel("Storage or Storage Deficit (mm)")
plt.show()
plt.savefig('FIG21.png', format='png', dpi=1600)

#%%



from math import isclose; import matplotlib.pyplot as plt
if not isclose(inset_ax._get_aspect_ratio(), ax3._get_aspect_ratio()):
    print("chosen inset x/ylim & width/height skew aspect w.r.t. main axes!")

ax[0].set_ylabel("Precipitation (mm/h)")

ax3.set_ylabel("Streamflow (mm/h)")


#%%
from math import isclose; import matplotlib.pyplot as plt

# set up main fig/axes
fig, main_ax = plt.subplots(); main_ax.set_box_aspect(0.5) 
inset_ax = main_ax.inset_axes(
   [0.05, 0.65, 0.3, 0.3],  # [x, y, width, height] w.r.t. axes
    xlim=[4, 5], ylim=[4, 5], # sets viewport & tells relation to main axes
    xticklabels=[], yticklabels=[]
)

# add plot content
for ax in main_ax, inset_ax:
    ax.plot([0, 9], [0, 9])  # first example line
    ax.plot([0, 9], [1, 8])  # second example line

# add zoom leaders
main_ax.indicate_inset_zoom(inset_ax, edgecolor="blue")

# careful! warn if aspect ratio of inset axes doesn't match main axes 
if not isclose(inset_ax._get_aspect_ratio(), main_ax._get_aspect_ratio()):
    print("chosen inset x/ylim & width/height skew aspect w.r.t. main axes!")
# %%
import matplotlib.pyplot as plt

from mpl_toolkits.axes_grid1.inset_locator import inset_axes

fig, (ax, ax2) = plt.subplots(1, 2, figsize=[5.5, 2.8])

# Create inset of width 1.3 inches and height 0.9 inches
# at the default upper right location
axins = inset_axes(ax, width=1.3, height=0.9)

# Create inset of width 30% and height 40% of the parent axes' bounding box
# at the lower left corner (loc=3)
axins2 = inset_axes(ax, width="30%", height="40%", loc=3)

# Create inset of mixed specifications in the second subplot;
# width is 30% of parent axes' bounding box and
# height is 1 inch at the upper left corner (loc=2)
axins3 = inset_axes(ax2, width="30%", height=1., loc=2)

# Create an inset in the lower right corner (loc=4) with borderpad=1, i.e.
# 10 points padding (as 10pt is the default fontsize) to the parent axes
axins4 = inset_axes(ax2, width="20%", height="20%", loc=4, borderpad=1)

# Turn ticklabels of insets off
for axi in [axins, axins2, axins3, axins4]:
    axi.tick_params(labelleft=False, labelbottom=False)

plt.show()

#%%
import numpy as np

from matplotlib import cbook
from matplotlib import pyplot as plt

fig, ax = plt.subplots()

# make data
Z = np.zeros((15, 15))  # 15x15 array
Z2 = np.zeros((150, 150))
ny, nx = Z.shape
Z2[30:30+ny, 30:30+nx] = Z
extent = (-3, 4, -4, 3)

ax.imshow(Z2, extent=extent, origin="lower")

# inset axes....
x1, x2, y1, y2 = -1.5, -0.9, -2.5, -1  # subregion of the original image
axins = ax.inset_axes(
    [0.5, 0.5, 0.47, 0.6],
    xlim=(x1, x2), ylim=(y1, y2), xticklabels=[], yticklabels=[])
axins.imshow(Z2, extent=extent, origin="lower")

ax.indicate_inset_zoom(axins, edgecolor="black")

plt.show()
# %%
import matplotlib.pyplot as plt
import numpy as np
from pandas import DataFrame
df = DataFrame(np.random.randn(5, 3), columns=['A', 'B', 'C'])

fig, ax = plt.subplots()
ax3 = ax.twinx()
rspine = ax3.spines['right']
rspine.set_position(('axes', 1.15))
ax3.set_frame_on(True)
ax3.patch.set_visible(False)
fig.subplots_adjust(right=0.7)

df.A.plot(ax=ax, style='b-')
# same ax as above since it's automatically added on the right
df.B.plot(ax=ax, style='r-', secondary_y=True)
df.C.plot(ax=ax3, style='g-')

# add legend --> take advantage of pandas providing us access
# to the line associated with the right part of the axis
ax3.legend([ax.get_lines()[0], ax.right_ax.get_lines()[0], ax3.get_lines()[0]],\
           ['A','B','C'], bbox_to_anchor=(1.5, 0.5))