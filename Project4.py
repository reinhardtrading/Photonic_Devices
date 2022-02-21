#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 25 13:19:46 2020

@author: bonhardt
"""


import numpy as np
import fd
import matplotlib.pyplot as plt
V=2.405
deltaepsilon =0.08;
na=np.sqrt(deltaepsilon);
a=1.5
wl=(2*np.pi*na*a)/V
print("cut-off wavelength is", wl)
r1=a;
na=V*wl/(2*np.pi*a)
ncl=np.sqrt(fd.sellmeier(wl))
nc=np.sqrt(na**2+ncl**2)
n1=nc;


#deltaepsilon=na**2
#print ("deltaepsilon is: ")
#print (deltaepsilon)
ncl=np.sqrt(fd.sellmeier(1.55))
nc=np.sqrt(na**2+ncl**2)
eps,wid=fd.epssif_tav(15,0.2,a,nc,ncl)
neff,hx,hy=fd.solve(1.55,eps,wid,1,nc,1,1)
neff1=neff
#print ("neff1: ")
#print (neff1)


ncl=np.sqrt(fd.sellmeier(1.551))
nc=np.sqrt(na**2+ncl**2)
eps,wid=fd.epssif_tav(15,0.2,a,nc,ncl)
neff,hx,hy=fd.solve(1.551,eps,wid,1,nc,1,1)
neff2=neff
#print ("neff2: ")
#print (neff2)


ncl=np.sqrt(fd.sellmeier(1.549))
nc=np.sqrt(na**2+ncl**2)
eps,wid=fd.epssif_tav(15,0.2,a,nc,ncl)
neff,hx,hy=fd.solve(1.549,eps,wid,1,nc,1,1)
neff0=neff
#print ("neff0: ")
#print (neff0)

D=(-1.55*(neff2+neff0-2*neff1))/(((0.001)**2)*2.99792e8)*1e12
print ("D is: ")
print (D)


#printing the slope
na=(2.405*1.55)/(2*np.pi*a)
deps=pow(na,2)
print("depsilon:",deps)
width=15
dx=0.2
a=1.5
wl=1.55
dwl=0.001
c=2.99792e8
n=0
neff=np.zeros(5)
wavelength = [wl-2*dwl,wl-dwl,wl,wl+dwl,wl+2*dwl]
for i in wavelength:
 ncl=np.sqrt(fd.sellmeier(i))
 nc=np.sqrt(pow(ncl,2)+deps)
 eps, wid = fd.epssif_tav(width, dx, a, nc, ncl)
 neff[n], hx, hy = fd.solve(i, eps, wid, 1, nc, 1, 1)
 n +=1
D=((-wl/c)*(neff[3]+neff[1]-2*neff[2])/pow(dwl,2))*1e12
print("D: ",D)
D1=((-(wl+dwl)/c)*(neff[4]+neff[2]-2*neff[3])/pow(dwl,2))*1e12
print("D1: ",D1)
D2=((-(wl-dwl)/c)*(neff[2]+neff[0]-2*neff[1])/pow(dwl,2))*1e12
print("D2: ",D2)
S=((D1-D2)/(2*dwl))*1e-3
print("Slope: ",S) 



#second part of the first excersise 
V=2.405
width=15;
dx=0.2
deltaepsilon =-0.03;
ncl=np.sqrt(fd.sellmeier(wl))
nco=np.sqrt((ncl**2)-deltaepsilon)
NA= np.sqrt((nco**2)-(ncl**2))
r2 = (V*wl)/(2*np.pi*NA)
print("R2 is", r2)


#na=V*wl/(2*np.pi*a)
ncl=np.sqrt(fd.sellmeier(wl))
nc=np.sqrt(NA**2+ncl**2)
#deltaepsilon=na**2
#print ("deltaepsilon is: ")
#print (deltaepsilon)
ncl=np.sqrt(fd.sellmeier(1.55))
nc=np.sqrt(NA**2+ncl**2)
n2=nc;

#eps,wid=fd.epssif_tav(15,0.2,r2,nc,ncl)
eps,wid=fd.epsdcf(width,dx,[r1,r2],[n1,n2,ncl])
neff,hx,hy=fd.solve(1.55,eps,wid,1,n2,1,1)
neff1=neff
#print ("neff1: ")
#print (neff1)


ncl=np.sqrt(fd.sellmeier(1.551))
nc=np.sqrt(NA**2+ncl**2)
n2=nc;

#eps,wid=fd.epssif_tav(15,0.2,r2,nc,ncl)
eps,wid=fd.epsdcf(width,dx,[r1,r2],[n1,n2,ncl])
neff,hx,hy=fd.solve(1.551,eps,wid,1,n2,1,1)
neff2=neff
#print ("neff2: ")
#print (neff2)



ncl=np.sqrt(fd.sellmeier(1.549))
nc=np.sqrt(NA**2+ncl**2)

n2=nc;
#eps,wid=fd.epssif_tav(15,0.2,r2,nc,ncl)
eps,wid=fd.epsdcf(width,dx,[r1,r2],[n1,n2,ncl])
neff,hx,hy=fd.solve(1.549,eps,wid,1,n2,1,1)
neff0=neff
#print ("neff0: ")
#print (neff0)
Dc=(-1.55*(neff2+neff0-2*neff1))/(((0.001)**2)*2.99792e8)*1e12
print ("Dispersion of the compensating fiber is: ")
print (Dc)











