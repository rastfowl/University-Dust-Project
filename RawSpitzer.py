# -*- coding: utf-8 -*-
"""
Created on Mon Oct 28 08:19:12 2019

@author: rastf
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Oct 21 07:26:54 2019

@author: rastf
"""

from astropy.io import fits
from astropy.wcs import WCS
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec
from matplotlib.colors import LogNorm
from matplotlib.patches import Circle
from astropy.visualization import simple_norm

#finding out what the coordinate system for each image is
wcs24 = WCS(fits.getheader('G54_Spitzer_24.fits'))
wcs70 = WCS(fits.getheader('G54_PACS_70.fits'))
wcs160 = WCS(fits.getheader('G54_PACS_160.fits'))
wcs250 = WCS(fits.getheader('G54_SPIRE_250.fits'))
wcs350 = WCS(fits.getheader('G54_SPIRE_350.fits'))
wcs500 = WCS(fits.getheader('G54_SPIRE_500.fits'))
wcs870 = WCS(fits.getheader('G54.1_FINAL_870.fits'))

#setting up the grid for plotting
f= plt.figure(figsize = (5,5))
gs1 = gridspec.GridSpec(1,1)
gs1.update(wspace=0.025, hspace=0)

#importing the data in
SPITZER24=fits.getdata('G54_Spitzer_24.fits')
PACS70=fits.getdata('G54_PACS_70.fits')
PACS160=fits.getdata('G54_PACS_160.fits')
SPIRE250=fits.getdata('G54_SPIRE_250.fits')
SPIRE350=fits.getdata('G54_SPIRE_350.fits')
SPIRE500=fits.getdata('G54_SPIRE_500.fits')
FINAL870=fits.getdata('G54.1_FINAL_870.fits')

#setting vmin and vmax
vmin=-0.03
vmax=1
cmap='cubehelix'
#defining some lists that will be needed for plotting
ax = []
lat = []
lon = []

# =============================================================================
# SPITZER 24 IMAGE
# =============================================================================
t=0
ax.append(plt.subplot(gs1[t], projection=wcs24)) #creating subplot and assigning correct wcs
img = ax[t].imshow(SPITZER24, cmap=cmap, origin = 'lower', norm=simple_norm(SPITZER24)) #plotting image with LogNorm
ax[t].text(0.5,0.92,'Spitzer MIPS 24$\mu$m', fontsize=13, family='serif', ha = 'left', va = 'center', transform = ax[t].transAxes, color='white') #adding text
c = Circle((292.6167, 18.8683), 0.025, edgecolor='white', facecolor='none',alpha=0.5,transform=ax[t].get_transform('fk5')) #defining circle
p = Circle((292.6292, 18.8667), 0.001, edgecolor='white', facecolor='none',alpha=0.5,transform=ax[t].get_transform('fk5')) #defining circle
ax[t].add_patch(c) #adding the circle
ax[t].add_patch(p)
#c2 = Circle((292.6167, 18.8683), 0.022, edgecolor='white', facecolor='none',alpha=0.5,transform=ax[t].get_transform('fk5')) #defining circle
#ax[t].add_patch(c2) #adding the circle2
lon.append(ax[t].coords[0]) #defining x axis
lat.append(ax[t].coords[1]) #defining y axis
lon[t].set_ticklabel_visible(True) #choosing to show or hide the x axis
lon[t].set_ticks_visible(True)
lat[t].set_ticklabel_visible(True) #choosing to show or hide the y axis
lat[t].set_ticks_visible(True)
ax[t].set_xlim((51,101)) #restricting field of view
ax[t].set_ylim((51,101)) #restricting field of view



