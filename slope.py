#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 25 13:09:47 2020

@author: bonhardt
"""


import fd
import numpy as np
#na=(2.405*1.1)/(2*np.pi*3.425)

deltaeps=-0.03
#deps=pow(na,2)
#print("depsilon:",deps)
width=15
dx=0.2
r1=1.5;
wl=1.55
na1=0.2828427
V=2.405;
deltaeps=-0.03
dwl=0.001
r2=3.425
c=2.99792e8
n=0

neff=np.zeros(5)
wavelength = [wl-2*dwl,wl-dwl,wl,wl+dwl,wl+2*dwl]
for i in wavelength:
 ncl=np.sqrt(fd.sellmeier(i))
 n1=np.sqrt(na1**2+ncl**2)
 nc=np.sqrt(pow(ncl,2)+deltaeps)
 n0=(nc+n1)/2
 #eps, wid = fd.epssif_tav(width, dx, a, nc, ncl)
 eps,wid=fd.epsdcf(width,dx,[r1,r2],[n1,nc,ncl])
 neff[n], hx, hy = fd.solve(i, eps, wid, 1, n0, 1, 1)
 n +=1
D=((-wl/c)*(neff[3]+neff[1]-2*neff[2])/pow(dwl,2))*1e12
print("D: ",D)
D1=((-(wl+dwl)/c)*(neff[4]+neff[2]-2*neff[3])/pow(dwl,2))*1e12
print("D1: ",D1)
D2=((-(wl-dwl)/c)*(neff[2]+neff[0]-2*neff[1])/pow(dwl,2))*1e12
print("D2: ",D2)
S=((D1-D2)/(2*dwl))*1e-3
print("Slope: ",S) 