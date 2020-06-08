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

#some extra stuff to enable automatic bounds for images
limits24 = WCS('G54_Spitzer_24.fits')
limits70 = WCS('G54_PACS_70.fits')
limits160 = WCS('G54_PACS_160.fits')
limits250 = WCS('G54_SPIRE_250.fits')
limits350 = WCS('G54_SPIRE_350.fits')
limits500 = WCS('G54_SPIRE_500.fits')
limits870 = WCS('G54.1_FINAL_870.fits')

#%%

#the bounds of each image are anchored by the pixel coordinates (50,50) and (100,100) in the PACS 70 image
baselftbound1=limits70.all_pix2world(49,50,0)[0]
baselftbound2=limits70.all_pix2world(51,50,0)[0]
basergtbound1=limits70.all_pix2world(99,50,0)[0]
basergtbound2=limits70.all_pix2world(101,50,0)[0]
basetopbound1=limits70.all_pix2world(50,99,0)[1]
basetopbound2=limits70.all_pix2world(50,101,0)[1]
basebotbound1=limits70.all_pix2world(50,49,0)[1]
basebotbound2=limits70.all_pix2world(50,51,0)[1]


for i in range(SPIRE250.shape[0]+1): #the limits for this loop at that of SPIRE250 as it had the largest dimensions
    
    if (baselftbound1>limits70.all_pix2world(i,50,0)[0])&(limits70.all_pix2world(i,50,0)[0]>baselftbound2):
        leftbound70=i
    if (basergtbound1>limits70.all_pix2world(i,50,0)[0])&(limits70.all_pix2world(i,50,0)[0]>basergtbound2):
        rightbound70=i
    if (basebotbound1<limits70.all_pix2world(50,i,0)[1])&(limits70.all_pix2world(50,i,0)[1]<basebotbound2):
        bottombound70=i
    if (basetopbound1<limits70.all_pix2world(50,i,0)[1])&(limits70.all_pix2world(50,i,0)[1]<basetopbound2):
        topbound70=i

#PACS 160 bounds
    if (baselftbound1>limits160.all_pix2world(i,50,0)[0])&(limits160.all_pix2world(i,50,0)[0]>baselftbound2):
        leftbound160=i
    if (basergtbound1>limits160.all_pix2world(i,50,0)[0])&(limits160.all_pix2world(i,50,0)[0]>basergtbound2):
        rightbound160=i
    if (basebotbound1<limits160.all_pix2world(50,i,0)[1])&(limits160.all_pix2world(50,i,0)[1]<basebotbound2):
        bottombound160=i
    if (basetopbound1<limits160.all_pix2world(50,i,0)[1])&(limits160.all_pix2world(50,i,0)[1]<basetopbound2):
        topbound160=i

#SPIRE 250 bounds
    if (baselftbound1>limits250.all_pix2world(i,50,0)[0])&(limits250.all_pix2world(i,50,0)[0]>baselftbound2):
        leftbound250=i
    if (basergtbound1>limits250.all_pix2world(i,50,0)[0])&(limits250.all_pix2world(i,50,0)[0]>basergtbound2):
        rightbound250=i
    if (basebotbound1<limits250.all_pix2world(50,i,0)[1])&(limits250.all_pix2world(50,i,0)[1]<basebotbound2):
        bottombound250=i
    if (basetopbound1<limits250.all_pix2world(50,i,0)[1])&(limits250.all_pix2world(50,i,0)[1]<basetopbound2):
        topbound250=i
        
#SPIRE 350 bounds
    if (baselftbound1>limits350.all_pix2world(i,50,0)[0])&(limits350.all_pix2world(i,50,0)[0]>baselftbound2):
        leftbound350=i
    if (basergtbound1>limits350.all_pix2world(i,50,0)[0])&(limits350.all_pix2world(i,50,0)[0]>basergtbound2):
        rightbound350=i
    if (basebotbound1<limits350.all_pix2world(50,i,0)[1])&(limits350.all_pix2world(50,i,0)[1]<basebotbound2):
        bottombound350=i
    if (basetopbound1<limits350.all_pix2world(50,i,0)[1])&(limits350.all_pix2world(50,i,0)[1]<basetopbound2):
        topbound350=i

#SPIRE 500 bounds
    if (baselftbound1>limits500.all_pix2world(i,50,0)[0])&(limits500.all_pix2world(i,50,0)[0]>baselftbound2):
        leftbound500=i
    if (basergtbound1>limits500.all_pix2world(i,50,0)[0])&(limits500.all_pix2world(i,50,0)[0]>basergtbound2):
        rightbound500=i
    if (basebotbound1<limits500.all_pix2world(50,i,0)[1])&(limits500.all_pix2world(50,i,0)[1]<basebotbound2):
        bottombound500=i
    if (basetopbound1<limits500.all_pix2world(50,i,0)[1])&(limits500.all_pix2world(50,i,0)[1]<basetopbound2):
        topbound500=i
        
#LABOCA 870 bounds
    if (baselftbound1>limits870.all_pix2world(i,50,0)[0])&(limits870.all_pix2world(i,50,0)[0]>baselftbound2):
        leftbound870=i
    if (basergtbound1>limits870.all_pix2world(i,50,0)[0])&(limits870.all_pix2world(i,50,0)[0]>basergtbound2):
        rightbound870=i
    if (basebotbound1<limits870.all_pix2world(50,i,0)[1])&(limits870.all_pix2world(50,i,0)[1]<basebotbound2):
        bottombound870=i
    if (basetopbound1<limits870.all_pix2world(50,i,0)[1])&(limits870.all_pix2world(50,i,0)[1]<basetopbound2):
        topbound870=i
#%%
    
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
ax[t].set_xlim(leftbound70,rightbound70) #restricting field of view
ax[t].set_ylim(bottombound70,topbound70) #restricting field of view



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
ax[t].set_xlim(leftbound160,rightbound160) #restricting field of view
ax[t].set_ylim(bottombound160,topbound160) #restricting field of view



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
ax[t].set_xlim(leftbound250,rightbound250) #restricting field of view
ax[t].set_ylim(bottombound250,topbound250) #restricting field of view


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
ax[t].set_xlim(leftbound350,rightbound350) #restricting field of view
ax[t].set_ylim(bottombound350,topbound350) #restricting field of view

# =============================================================================
# SPIRE 500 IMAGE - 43x43
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
ax[t].set_xlim(leftbound500,rightbound500) #restricting field of view
ax[t].set_ylim(bottombound500,topbound500) #restricting field of view

# =============================================================================
# LABOCA 870 IMAGE - 75x75
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
ax[t].set_xlim(leftbound870,rightbound870) #restricting field of view
ax[t].set_ylim(bottombound870,topbound870) #restricting field of view


