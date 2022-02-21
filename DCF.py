#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 26 12:17:51 2020

@author: bonhardt
"""


import numpy as np
import fd
#import matplotlib.pyplot as plt

#deltaepssilon=0.08
width=15
dx=0.2
r1=1.5;
wl=1.55
na1=0.2828427
V=2.405;
deltaeps=-0.03
dl=0.001


ncl=np.sqrt(fd.sellmeier(wl))
nco=np.sqrt(ncl**2+deltaeps)
#NA= np.sqrt((nco**2)-(ncl**2))
NA=np.sqrt(-deltaeps)
r2 = (V*wl)/(2*np.pi*NA)
print("R2 is", r2)

#neff of wl
ncl1=np.sqrt(fd.sellmeier(wl));
n2=np.sqrt(ncl1**2+deltaeps)
n1=np.sqrt(na1**2+ncl1**2)
nc=(n1+n2)/2
eps,wid=fd.epsdcf(width,dx,[r1,r2],[n1,n2,ncl1])
neff,hx,hy=fd.solve(wl,eps,wid,1,nc,1,1)
neff0=neff

#neff of wl-dl
ncl2=np.sqrt(fd.sellmeier(wl-dl));
n2=np.sqrt(ncl2**2+deltaeps)
n1=np.sqrt(na1**2+ncl2**2)
nc=(n1+n2)/2

eps,wid=fd.epsdcf(width,dx,[r1,r2],[n1,n2,ncl2])
neff,hx,hy=fd.solve(wl-dl,eps,wid,1,nc,1,1)
neff1=neff

#neff of wl+dl
ncl3=np.sqrt(fd.sellmeier(wl+dl));
n2=np.sqrt(ncl3**2+deltaeps)
n1=np.sqrt(na1**2+ncl3**2)
nc=(n1+n2)/2

eps,wid=fd.epsdcf(width,dx,[r1,r2],[n1,n2,ncl3])
neff,hx,hy=fd.solve(wl+dl,eps,wid,1,nc,1,1)
neff2=neff

Dc0=(-1.55*(neff1+neff2-2*neff0))/(((0.001)**2)*2.99792e8)*1e12
print ("Dispersion of the compensating fiber at 1.55um is: ")
print (Dc0)


#DCF at 1.530 


#neff of wl

wl =1.53;
ncl1=np.sqrt(fd.sellmeier(wl));
n2=np.sqrt(ncl1**2+deltaeps)
n1=np.sqrt(na1**2+ncl1**2)
nc=(n1+n2)/2
eps,wid=fd.epsdcf(width,dx,[r1,r2],[n1,n2,ncl1])
neff,hx,hy=fd.solve(wl,eps,wid,1,nc,1,1)
neff0=neff

#neff of wl-dl
ncl2=np.sqrt(fd.sellmeier(wl-dl));
n2=np.sqrt(ncl2**2+deltaeps)
n1=np.sqrt(na1**2+ncl2**2)
nc=(n1+n2)/2

eps,wid=fd.epsdcf(width,dx,[r1,r2],[n1,n2,ncl2])
neff,hx,hy=fd.solve(wl-dl,eps,wid,1,nc,1,1)
neff1=neff

#neff of wl+dl
ncl3=np.sqrt(fd.sellmeier(wl+dl));
n2=np.sqrt(ncl3**2+deltaeps)
n1=np.sqrt(na1**2+ncl3**2)
nc=(n1+n2)/2

eps,wid=fd.epsdcf(width,dx,[r1,r2],[n1,n2,ncl3])
neff,hx,hy=fd.solve(wl+dl,eps,wid,1,nc,1,1)
neff2=neff

Dc1=(-wl*(neff1+neff2-2*neff0))/(((0.001)**2)*2.99792e8)*1e12
print ("Dispersion of the compensating fiber at 1.53um is: ")
print (Dc1)



wl=1.565
#neff of wl
ncl1=np.sqrt(fd.sellmeier(wl));
n2=np.sqrt(ncl1**2+deltaeps)
n1=np.sqrt(na1**2+ncl1**2)
nc=(n1+n2)/2
eps,wid=fd.epsdcf(width,dx,[r1,r2],[n1,n2,ncl1])
neff,hx,hy=fd.solve(wl,eps,wid,1,nc,1,1)
neff0=neff

#neff of wl-dl
ncl2=np.sqrt(fd.sellmeier(wl-dl));
n2=np.sqrt(ncl2**2+deltaeps)
n1=np.sqrt(na1**2+ncl2**2)
nc=(n1+n2)/2

eps,wid=fd.epsdcf(width,dx,[r1,r2],[n1,n2,ncl2])
neff,hx,hy=fd.solve(wl-dl,eps,wid,1,nc,1,1)
neff1=neff

#neff of wl+dl
ncl3=np.sqrt(fd.sellmeier(wl+dl));
n2=np.sqrt(ncl3**2+deltaeps)
n1=np.sqrt(na1**2+ncl3**2)
nc=(n1+n2)/2

eps,wid=fd.epsdcf(width,dx,[r1,r2],[n1,n2,ncl3])
neff,hx,hy=fd.solve(wl+dl,eps,wid,1,nc,1,1)
neff2=neff

Dc2=(-wl*(neff1+neff2-2*neff0))/(((0.001)**2)*2.99792e8)*1e12
print ("Dispersion of the compensating fiber at 1.565um is: ")
print (Dc2)

DCF = Dc0;          #DCF for 1.55um
DCF1= Dc1;          #DCF for 1.53um
DCF2=Dc2;           #DCF for 1.565um
LTF=100
DTF=17.48
LDC=-(LTF*DTF)/DCF
print("The length of compensating fibre is", LDC);

Tot1= (LTF*DTF)+(LDC*DCF1)
Tot2= (LTF*DTF)+(LDC*DCF2)

print("The total dispersion at 1.53um is", Tot1)
print("The total dispersion at 1.565um is", Tot2)





