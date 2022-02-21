#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 17 20:40:24 2020

@author: bonhardt
"""


#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 15 11:46:26 2020

@author: bonhardt
"""


import numpy as np
import matplotlib.pyplot as plt
import fd


a= 4.5;
V = 2.405;
wl=1.3;
width =15;
dx=0.2;
NA=(V*wl)/(2*np.pi*a);
deltaepsilon = NA**2;
print ("delta epsilon is", deltaepsilon);



ncl=np.sqrt(fd.sellmeier(1.3))
nco =np.sqrt(NA**2+ncl**2);
#print ("ncore:", nco)
#print("ncl", ncl)
eps,wid=fd.epssif_tav(width,dx,a,nco,ncl,quarter=0)
neff0,hx,hy=fd.solve(1.3,eps,wid,1,nco,1,1)
neff_0= neff0;


ncl1=np.sqrt(fd.sellmeier(1.8))
nco1 =np.sqrt(NA**2+ncl1**2);
#print ("ncore1:", nco1)
#print("ncl1", ncl1)
eps,wid=fd.epssif_tav(width,dx,a,nco1,ncl1,quarter=0)
neff2,hx,hy=fd.solve(1.8,eps,wid,1,nco1,1,1)
neff_2= neff2;

ncl2=np.sqrt(fd.sellmeier(1.55))
nco2 =np.sqrt(NA**2+ncl2**2);
#print ("ncore2:", nco1)
#print("ncl2", ncl1)
eps,wid=fd.epssif_tav(width,dx,a,nco2,ncl2,quarter=0)
neff1,hx,hy=fd.solve(1.55,eps,wid,1,nco2,1,1)
neff_1= 2*neff1;


second_derivative= ((neff_0+neff_2)-neff_1)/(0.25*0.25);

D= -(1.55*second_derivative)/(3.0*1e8);
Dispersion = D*1e12;
print(Dispersion)
