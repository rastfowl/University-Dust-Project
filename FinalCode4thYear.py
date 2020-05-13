# -*- coding: utf-8 -*-
"""
Created on Fri Feb 21 12:12:37 2020

@author: rastfowl
"""
#preamble
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.visualization import simple_norm
from photutils import SkyCircularAperture, aperture_photometry, CircularAnnulus
from matplotlib.colors import LogNorm, Normalize
from matplotlib.patches import Circle
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import astropy.units as u
from tabulate import tabulate

#load in data
datapath='G54_tdenscube.fits'
headpath='G54_sigtdenscube.fits'

#this is for latex output
latex=np.zeros((4,12,5))
MASSES=np.zeros(4)
ERROR=np.zeros(4)

#defining variables and setting up temp array
betas=np.array([1.5,1.75,2.00,2.25])
tempno=12
lowertemp=20.
highertemp=90.
temp = np.logspace(np.log10(lowertemp),np.log10(highertemp),tempno)
smallradius=0.025
largeradius=0.03
H2mass=(2.72*1.67e-27)
anglesubtendedby1pix=2e-5
distanceincm=1.851e22
solmass=1.989e30
rightascension=292.6167
declination=18.8683

#loading in data 
dataA = fits.getdata(datapath)[0,:,:,:]
dataB = fits.getdata(datapath)[1,:,:,:]
dataC = fits.getdata(datapath)[2,:,:,:]
dataD = fits.getdata(datapath)[3,:,:,:]
bigdata=dataA,dataB,dataC,dataD

#converting data from Gas Mass (10^20H2/cm2) to Dust mass (M_(odot)/pix)
H2percm2=1e20
kgpercm2=H2percm2*H2mass
cmin1pix=anglesubtendedby1pix*distanceincm
cm2in1pix=cmin1pix**2
dr=anglesubtendedby1pix*206265
kgperpix=kgpercm2*cm2in1pix
solsperpix=kgperpix/solmass
dustsolsperpix=solsperpix/100

#converting data to solar mass units
dataAc=dataA*dustsolsperpix
dataBc=dataB*dustsolsperpix
dataCc=dataC*dustsolsperpix
dataDc=dataD*dustsolsperpix
bigdatac=dataAc,dataBc,dataCc,dataDc

#creating blank arrays for filling up later
M_S=np.zeros((4,tempno)) #mass in small ap
M_B=np.zeros((4,tempno)) #mass in large ap
M_bg=np.zeros((4,tempno)) #mass in annulus between
M_subt=np.zeros((4,tempno)) #mass in annulus between
rho_bg=np.zeros((4,tempno)) #background mass density
error=np.zeros((4,tempno)) #error in mass of small ap

#number of pixels in each aperture
A_1pix=4**2
A_S=np.pi*((smallradius*3600)**2)
A_B=np.pi*((largeradius*3600)**2)
N_S=A_S/A_1pix
N_B=A_B/A_1pix
N_bg=N_B-N_S

#setting up the loop that will be run through four times for each beta
for i in range(0,4):
    sample=bigdata[i]
    
    wcs = WCS(fits.getheader(headpath)).sub(axes=2)
    norm = simple_norm(sample, 'sqrt', max_percent=88)
    figure = plt.figure(figsize=(9,7))
    grid = gridspec.GridSpec(3,4, wspace=0, hspace=0)
    ax = []
    lat = []
    lon = []
        
    #the loop that runs through each temperature map
    for t in range(tempno):
        
        #plotting
        ax.append(plt.subplot(grid[t], projection=wcs)) #defining the subplot for this iteration
        img=ax[t].imshow(sample[t,:,:],cmap='cubehelix',origin='lower',norm=norm) #plotting data
        ax[t].text(0.66,0.1,'{:0.1f} K'.format(temp[t]), \
        transform=ax[t].transAxes, color='black', bbox=dict(facecolor='white')) #adding temperature caption
        smallcirc=Circle((rightascension,declination),smallradius,edgecolor='cyan',facecolor='none',alpha=0.5, \
                         transform=ax[t].get_transform('fk5')) #small circle
        largecirc=Circle((rightascension,declination),largeradius,edgecolor='cyan',facecolor='none',alpha=0.5, \
                         transform=ax[t].get_transform('fk5')) #large circle
        ax[t].add_patch(smallcirc)
        ax[t].add_patch(largecirc)
        #plt.suptitle('Temperature grid for \u03B2 = %.3f' % betas[i]) #optional title
        
        #plot formatting
        lon.append(ax[t].coords[0])
        lat.append(ax[t].coords[1])
        ax[t].set_xlim((25, 95))
        ax[t].set_ylim((15, 85))
        
        #hiding axes labels
        lon[t].set_ticklabel_visible(False)
        lat[t].set_ticklabel_visible(False)
        lon[t].set_ticks_visible(False)
        lat[t].set_ticks_visible(False)
        
        #this makes it so that only the outer plots have labels
        if t==0.00:
            lon[t].set_ticklabel_visible(False)
            lat[t].set_ticklabel_visible(True)
            lat[t].set_ticks_visible(True)
            lat[t].set_ticks_position('l')
        if t==4.00:
            lon[t].set_ticklabel_visible(False)
            lat[t].set_ticklabel_visible(True)
            lat[t].set_ticks_visible(True)
            lat[t].set_ticks_position('l')
        if t==8.00:
            lon[t].set_ticklabel_visible(True)
            lat[t].set_ticklabel_visible(True)
            lat[t].set_ticks_visible(True)
            lat[t].set_ticks_position('l')
            lon[t].set_ticks_visible(True)
            lon[t].set_ticks_position('b')
            lon[t].set_ticklabel(exclude_overlapping=True)
        if t==9.00:
            lon[t].set_ticklabel_visible(True)
            lat[t].set_ticklabel_visible(False)
            lon[t].set_ticks_visible(True)
            lon[t].set_ticks_position('b')
            lon[t].set_ticklabel(exclude_overlapping=True)
        if t==10.00:
            lon[t].set_ticklabel_visible(True)
            lat[t].set_ticklabel_visible(False)
            lon[t].set_ticks_visible(True)
            lon[t].set_ticks_position('b')
            lon[t].set_ticklabel(exclude_overlapping=True)
        if t==11.00:
            lon[t].set_ticklabel_visible(True)
            lat[t].set_ticklabel_visible(False)
            lon[t].set_ticks_visible(True)
            lon[t].set_ticks_position('b')
            lon[t].set_ticklabel(exclude_overlapping=True)
        
        #aperture photometry
        converted=sample[t,:,:]*dustsolsperpix
        positions=SkyCoord(ra=[rightascension],dec=[declination],unit='deg')
        smallaperture=SkyCircularAperture(positions, smallradius*u.deg)
        largeaperture=SkyCircularAperture(positions, largeradius*u.deg)
        smallphot_table=aperture_photometry(converted,smallaperture, wcs=wcs)
        largephot_table=aperture_photometry(converted,largeaperture, wcs=wcs)
        M_S[i,t]=(smallphot_table['aperture_sum']) #mass in small aperture
        M_B[i,t]=(largephot_table['aperture_sum']) #mass in large aperture
        
        #mass in background annulus
        M_bg[i,t]=M_B[i,t]-M_S[i,t]
        #background mass density
        rho_bg[i,t]=M_bg[i,t]/N_bg
        #background subtracted mass in small aperture
        M_subt[i,t]=M_S[i,t]-rho_bg[i,t]*N_S

        #errors
        snr_xposition=(smallphot_table['xcenter']).value[0]
        snr_yposition=(smallphot_table['ycenter']).value[0]
        annupositions=[snr_xposition,snr_yposition]
        annulus_aperture=CircularAnnulus(annupositions, r_in=(smallradius*3600)/4, r_out=(largeradius*3600)/4)
        annulus_masks=annulus_aperture.to_mask(method='center')
        annulus_data=annulus_masks.multiply(sample[t,:,:])
        mask=annulus_masks.data
        annulus_data_1d=annulus_data[mask>0]
        annulus_data_1d.shape
        stdev=np.std(annulus_data_1d*dustsolsperpix)
        ksi=0.1*M_S[i,t]
        error[i,t]=np.sqrt(0+((N_S**2)*(stdev**2/N_bg))+ksi**2)
    
    #colour bar
    cbar_ax=figure.add_axes([0.91, 0.11, 0.03, 0.765]) #Adding a colour bar
    figure.colorbar(img, cax=cbar_ax, orientation='vertical')
    
    #axes labels for plot
    figure.text(0.50, 0.04, 'Right Ascension $(\u03B1_{2000})$', ha='center')
    figure.text(0.04, 0.50, 'Declination $(\u03B4_{2000})$', va='center', rotation='vertical')
    figure.text(0.91, 0.04, '$M_{\odot} pix^{-1}$')
    
    #this prints out the results in the console
    headers=["Temp"]
    table=np.zeros((12,11))
    table[:,0]=temp
    headers.append("M_S")
    table[:,1]=M_S[i,:]
    headers.append("M_B")
    table[:,2]=M_B[i,:]
    headers.append("M_bg")
    table[:,3]=M_B[i,:]-M_S[i,:]
    headers.append("M/pix (S)")
    table[:,4]=M_S[i,:]/N_S
    headers.append("M/pix(bg)")
    table[:,5]=rho_bg[i,:]
    headers.append("M_S")
    table[:,6]=M_S[i,:]
    headers.append("M_bg")
    table[:,7]=rho_bg[i,:]*N_S
    headers.append("M_subt")
    table[:,8]=M_S[i,:]-rho_bg[i,:]*N_S
    headers.append("Error")
    table[:,9]=error[i,:]
    headers.append("% Error")
    table[:,10]=error[i,:]*100/(M_S[i,:]-rho_bg[i,:]*N_S)
    print('\n\n\n \033[4m This is the information for \u03B2 = %.2f \033[0m' % betas[i])
    print(tabulate(table, headers=headers, tablefmt="github"))
    
    #this records the total mass at each temperature
    MASSES[i]=np.sum(M_subt[i,:])
    ERRORsqrd=np.sum(error[i,:]**2)
    ERROR[i]=np.sqrt(ERRORsqrd)
    perc=ERROR[i]*100/MASSES[i]
    print("For \u03B2 = %.2f the mass is %.5f +/- %.5f" % (betas[i],MASSES[i],ERROR[i]))
    
    #this bit exports the array data to an external .txt file that I can import into my report
    latex[i,:,0]=temp
    latex[i,:,1]=M_S[i,:]
    latex[i,:,2]=rho_bg[i,:]
    latex[i,:,3]=M_S[i,:]-rho_bg[i,:]*N_S
    latex[i,:,4]=error[i,:]
    if i==0:
        string1='beta1dot5.txt'
    if i==1:
        string1='beta1dot75.txt'
    if i==2:
        string1='beta2dot00.txt'
    if i==3:
        string1='beta2dot25.txt'
    np.savetxt(string1,latex[i,:,:],delimiter='&',newline='\\\\', fmt='%.6e')