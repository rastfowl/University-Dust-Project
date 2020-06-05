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
f= plt.figure(figsize = (9,6))
gs1 = gridspec.GridSpec(2,3)
gs1.update(wspace=0.0, hspace=0.0)

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
# =============================================================================
# t=0
# ax.append(plt.subplot(gs1[t], projection=wcs24)) #creating subplot and assigning correct wcs
# img = ax[t].imshow(SPITZER24, cmap=cmap, origin = 'lower', norm=simple_norm(SPITZER24)) #plotting image with LogNorm
# ax[t].text(0.40,0.85,'Spitzer 24um', ha = 'left', va = 'center', transform = ax[t].transAxes, color='white') #adding text
# c = Circle((292.6167, 18.8683), 0.025, edgecolor='white', facecolor='none',alpha=0.5,transform=ax[t].get_transform('fk5')) #defining circle
# p = Circle((292.6292, 18.8667), 0.001, edgecolor='white', facecolor='none',alpha=0.5,transform=ax[t].get_transform('fk5')) #defining circle
# ax[t].add_patch(c) #adding the circle
# ax[t].add_patch(p)
# #c2 = Circle((292.6167, 18.8683), 0.022, edgecolor='white', facecolor='none',alpha=0.5,transform=ax[t].get_transform('fk5')) #defining circle
# #ax[t].add_patch(c2) #adding the circle2
# lon.append(ax[t].coords[0]) #defining x axis
# lat.append(ax[t].coords[1]) #defining y axis
# lon[t].set_ticklabel_visible(False) #choosing to show or hide the x axis
# lon[t].set_ticks_visible(False)
# lat[t].set_ticklabel_visible(False) #choosing to show or hide the y axis
# lat[t].set_ticks_visible(False)
# ax[t].set_xlim((30,110)) #restricting field of view
# ax[t].set_ylim((30,110)) #restricting field of view
# =============================================================================



# =============================================================================
# PACS 70 IMAGE - 50x50
# =============================================================================
t=0
ax.append(plt.subplot(gs1[t], projection=wcs70)) #creating subplot and assigning correct wcs
img = ax[t].imshow(PACS70, cmap=cmap, origin = 'lower',vmin = 0.01, vmax = 1.0, norm = simple_norm(PACS70, max_percent=90, stretch='sqrt')) #plotting image with LogNorm
ax[t].text(0.23,0.91,'Herschel PACS 70$\mu$m', fontsize=11, family='serif', ha = 'left', va = 'center', transform = ax[t].transAxes, color='white') #adding text
c = Circle((292.6167, 18.8683), 0.025, edgecolor='white', facecolor='none',alpha=0.5,transform=ax[t].get_transform('fk5')) #defining circle
p = Circle((292.6292, 18.8667), 0.001, edgecolor='white', facecolor='none',alpha=0.5,transform=ax[t].get_transform('fk5')) #defining circle
ax[t].add_patch(c) #adding the circle
ax[t].add_patch(p)
#c2 = Circle((292.6167, 18.8683), 0.022, edgecolor='white', facecolor='none',alpha=0.5,transform=ax[t].get_transform('fk5')) #defining circle
#ax[t].add_patch(c2) #adding the circle2
lon.append(ax[t].coords[0]) #defining x axis
lat.append(ax[t].coords[1]) #defining y axis
lon[t].set_ticklabel_visible(False) #choosing to show or hide the x axis
lon[t].set_ticks_visible(False)
lat[t].set_ticklabel_visible(True) #choosing to show or hide the y axis
lat[t].set_ticks_visible(True)
ax[t].set_xlim((50, 100)) #restricting field of view
ax[t].set_ylim((50, 100)) #restricting field of view



# =============================================================================
# PACS 160 IMAGE - 50x50
# =============================================================================
t=1
ax.append(plt.subplot(gs1[t], projection=wcs160)) #creating subplot and assigning correct wcs
img = ax[t].imshow(PACS160, cmap=cmap, origin = 'lower',vmin = 0.01, vmax = 1.0, norm = simple_norm(PACS160, max_percent=90, stretch='sqrt'))
ax[t].text(0.185,0.91,'Herschel PACS 160$\mu$m', fontsize=11, family='serif', ha = 'left', va = 'center', transform = ax[t].transAxes, color='white') #adding text
c = Circle((292.6167, 18.8683), 0.025, edgecolor='white', facecolor='none',alpha=0.5,transform=ax[t].get_transform('fk5')) #defining circle
p = Circle((292.6292, 18.8667), 0.001, edgecolor='white', facecolor='none',alpha=0.5,transform=ax[t].get_transform('fk5')) #defining circle
ax[t].add_patch(c) #adding the circle
ax[t].add_patch(p)
#c2 = Circle((292.6167, 18.8683), 0.022, edgecolor='white', facecolor='none',alpha=0.5,transform=ax[t].get_transform('fk5')) #defining circle
#ax[t].add_patch(c2) #adding the circle2
lon.append(ax[t].coords[0]) #defining x axis
lat.append(ax[t].coords[1]) #defining y axis
lon[t].set_ticklabel_visible(False) #choosing to show or hide the x axis
lon[t].set_ticks_visible(False)
lat[t].set_ticklabel_visible(False) #choosing to show or hide the y axis
lat[t].set_ticks_visible(False)
ax[t].set_xlim((50, 100)) #restricting field of view
ax[t].set_ylim((50, 100)) #restricting field of view



#repeat this for other images, but change LogNorm to something else to get the SPIRE images to plot. remember to put t=2 and t=3 etc.


# =============================================================================
# SPIRE 250 IMAGE - 100x100
# =============================================================================
t=2
ax.append(plt.subplot(gs1[t], projection=wcs250)) #creating subplot and assigning correct wcs
img = ax[t].imshow(SPIRE250, cmap=cmap, origin = 'lower', norm=simple_norm(SPIRE250, max_percent=90, stretch='power', power=2)) #plotting image with LogNorm
ax[t].text(0.155,0.91,'Herschel SPIRE 250$\mu$m', fontsize=11, family='serif', ha = 'left', va = 'center', transform = ax[t].transAxes, color='white') #adding text
c = Circle((292.6167, 18.8683), 0.025, edgecolor='white', facecolor='none',alpha=0.5,transform=ax[t].get_transform('fk5')) #defining circle
p = Circle((292.6292, 18.8667), 0.001, edgecolor='white', facecolor='none',alpha=0.5,transform=ax[t].get_transform('fk5')) #defining circle
ax[t].add_patch(c) #adding the circle
ax[t].add_patch(p)
#c2 = Circle((292.6167, 18.8683), 0.022, edgecolor='white', facecolor='none',alpha=0.5,transform=ax[t].get_transform('fk5')) #defining circle
#ax[t].add_patch(c2) #adding the circle2
lon.append(ax[t].coords[0]) #defining x axis
lat.append(ax[t].coords[1]) #defining y axis
lon[t].set_ticklabel_visible(False) #choosing to show or hide the x axis
lon[t].set_ticks_visible(False)
lat[t].set_ticklabel_visible(False) #choosing to show or hide the y axis
lat[t].set_ticks_visible(False)
ax[t].set_xlim((2790,2890)) #restricting field of view
ax[t].set_ylim((1076,1176)) #restricting field of view


# =============================================================================
# SPIRE 350 IMAGE - 60x60
# =============================================================================
t=3
ax.append(plt.subplot(gs1[t], projection=wcs350)) #creating subplot and assigning correct wcs
img = ax[t].imshow(SPIRE350, cmap=cmap, origin = 'lower', norm=simple_norm(SPIRE350, max_percent=90, stretch='power', power=2)) #plotting image with LogNorm
ax[t].text(0.155,0.91,'Herschel SPIRE 350$\mu$m', fontsize=11, family='serif', ha = 'left', va = 'center', transform = ax[t].transAxes, color='white') #adding text
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
ax[t].set_xlim((1677,1737)) #restricting field of view
ax[t].set_ylim((652,712)) #restricting field of view

# =============================================================================
# SPIRE 500 IMAGE - 44x44
# =============================================================================
t=4
ax.append(plt.subplot(gs1[t], projection=wcs500)) #creating subplot and assigning correct wcs
img = ax[t].imshow(SPIRE500, cmap=cmap, origin = 'lower', norm=simple_norm(SPIRE500,  max_percent=90, stretch='power', power=2)) #plotting image with LogNorm
ax[t].text(0.155,0.91,'Herschel SPIRE 500$\mu$m', fontsize=11, family='serif', ha = 'left', va = 'center', transform = ax[t].transAxes, color='white') #adding text
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
lat[t].set_ticklabel_visible(False) #choosing to show or hide the y axis
lat[t].set_ticks_visible(False)
ax[t].set_xlim((1198,1240)) #restricting field of view
ax[t].set_ylim((465,507)) #restricting field of view

# =============================================================================
# LABOCA 870 IMAGE
# =============================================================================
t=5
ax.append(plt.subplot(gs1[t], projection=wcs870)) #creating subplot and assigning correct wcs
img = ax[t].imshow(FINAL870, cmap=cmap, origin = 'lower', norm=simple_norm(FINAL870, min_percent=20, max_percent=90, stretch='power', power=2)) #plotting image with LogNorm
ax[t].text(0.18,0.91,'APEX LABOCA 870$\mu$m', fontsize=11, family='serif', ha = 'left', va = 'center', transform = ax[t].transAxes, color='white') #adding text
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
lat[t].set_ticklabel_visible(False) #choosing to show or hide the y axis
lat[t].set_ticks_visible(False)
ax[t].set_xlim(((128-36),(128+36))) #restricting field of view
ax[t].set_ylim(((136-36),(136+36))) #restricting field of view


