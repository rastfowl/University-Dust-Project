# -*- coding: utf-8 -*-
"""
Created on Wed May 13 12:06:05 2020

@author: rastf
"""
import matplotlib.pyplot as plt
import numpy as np

#%%
plt.figure()
    
#this cell plots the mass vs beta
p=np.polyfit(betas,MASSES,3)
z=np.poly1d(p)
x=np.linspace(1.5,2.25,100)

plt.plot(x,z(x), 'y--', label="Prediction")
plt.plot(betas,MASSES, '.', markersize=10, label="Data")
plt.errorbar(betas,MASSES, yerr=ERROR, fmt='none', color='orange', capsize=6)
plt.legend(loc="best")
plt.ylabel("$M_{G54.1}\,\,(M_{\odot})$")
plt.xlabel("Dust opacity exponent \u03B2")
plt.grid()

plt.savefig('massbetas2.png')

index=np.where(z(x)==max(z(x)))

maxweight=z(x)[index]
correspondingbeta=x[index]

#%%
plt.figure()
    
#this cell plots the 
p=np.polyfit(betas,masssums,3)
z=np.poly1d(p)
x=np.linspace(1.5,2.25,100)


plt.plot(x,z(x), 'y--', label="Prediction")
plt.plot(betas,masssums, '.', markersize=10, label="Data")
plt.errorbar(betas,masssums, yerr=erra, fmt='none', color='orange', capsize=6)
plt.legend(loc="best")
plt.ylabel("$M_{ISM}\,\,(M_{\odot})$")
plt.xlabel("Dust opacity exponent \u03B2")
plt.grid()

plt.savefig('rship2.png')


#%%

plt.figure(figsize=(10,7))

plt.subplot(221)
plt.plot(temp,M_S[0,:], 'r--', label='$M_{Source}$')
plt.plot(temp,M_ism[0,:], 'r-.', label='$M_{ISM}$') #this one
plt.plot(temp,M_subt[0,:], 'r-', label='$M_{G54.1}$')
plt.title('\u03B2=1.50')
plt.xlabel('Temperature $(K)$')
plt.ylabel('Mass $ (M_{\odot})$')
plt.legend(loc='best', fontsize=12)
plt.grid()

plt.subplot(222)
plt.plot(temp,M_S[1,:], 'y--', label='$M_{Source}$')
plt.plot(temp,M_ism[1,:], 'y-.', label='$M_{ISM}$') #this one
plt.plot(temp,M_subt[1,:], 'y-', label='$M_{G54.1}$')
plt.title('\u03B2=1.75')
plt.xlabel('Temperature $(K)$')
plt.ylabel('Mass $(M_{\odot})$')
plt.legend(loc='best', fontsize=12)
plt.grid()

plt.subplot(223)
plt.plot(temp,M_S[2,:], 'b--', label='$M_{Source}$')
plt.plot(temp,M_ism[2,:], 'b-.', label='$M_{ISM}$') #this one
plt.plot(temp,M_subt[2,:], 'b-', label='$M_{G54.1}$')
plt.title('\u03B2=2.00')
plt.xlabel('Temperature $(K)$')
plt.ylabel('Mass $(M_{\odot})$')
plt.legend(loc='best', fontsize=12)
plt.grid()

plt.subplot(224)
plt.plot(temp,M_S[3,:], 'g--', label='$M_{Source}$')
plt.plot(temp,M_ism[3,:], 'g-.', label='$M_{ISM}$') #this one
plt.plot(temp,M_subt[3,:], 'g-', label='$M_{G54.1}$')
plt.title('\u03B2=2.25')
plt.xlabel('Temperature $(K)$')
plt.ylabel('Mass $(M_{\odot})$')
plt.legend(loc='best', fontsize=12)
plt.grid()


plt.tight_layout()

plt.savefig('finalplot.png')


#%%
plt.figure(figsize=(10,7))

plt.subplot(221)
plt.plot(temp,M_S[0,:], 'r-', label='$M_{Source}$')

plt.title('\u03B2=1.50')
plt.xlabel('Temperature $(K)$')
plt.ylabel('Mass $ (M_{\odot})$')
plt.legend(loc='best')
plt.grid()

plt.subplot(222)
plt.plot(temp,M_S[1,:], 'y-', label='$M_{Source}$')

plt.title('\u03B2=1.75')
plt.xlabel('Temperature $(K)$')
plt.ylabel('Mass $(M_{\odot})$')
plt.legend(loc='best')
plt.grid()

plt.subplot(223)
plt.plot(temp,M_S[2,:], 'b-', label='$M_{Source}$')

plt.title('\u03B2=2.00')
plt.xlabel('Temperature $(K)$')
plt.ylabel('Mass $(M_{\odot})$')
plt.legend(loc='best')
plt.grid()

plt.subplot(224)
plt.plot(temp,M_S[3,:], 'g-', label='$M_{Source}$')

plt.title('\u03B2=2.25')
plt.xlabel('Temperature $(K)$')
plt.ylabel('Mass $(M_{\odot})$')
plt.legend(loc='best')
plt.grid()

plt.tight_layout()