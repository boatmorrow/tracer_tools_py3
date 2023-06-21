#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 15 14:31:39 2023

Test He_tools 39Ar production

@author: wpgardner
"""

import numpy as np
from tracer_tools.noble_gas_tools import lapse_rate
import pdb
import matplotlib.pyplot as plt
import tracer_tools.He_tools as he



#rock type
UC = he.rock_type()
UC.defined_rock_type()

#benchmark against JFM - test single depth
He = 0
Hd = 80
depth = 0.5#mwe
depth = 0.5*100 #1000k/m3 * 1/10 kg/m2 -> g/cm2
UC.APR.depth=depth
UC.APR.elev=He
UC.APR.incl=Hd
Pn = UC.calcPn_ev()
print('P_n = %1.2f n/g_rock/yr'%UC.Pn_ev)

UC.calcAr_prod_rate_whole_rock()
print('Total 39Ar production rate is %1.3g,\
 from hadronic flux %1.3g,\
 from muon adsorb. is %1.3g'\
          %(UC.APR.Ar39_value,UC.APR.P39Ar_ev_n,UC.APR.P39Ar_mu))

#heidelberg depth profile
He = 100 #elevation
Hd = 65 #magnetic inclination
depth_r = np.linspace(0,50) #m
depth_rw = depth_r*250 #assuming 2500kg/m3 rock density
UC.APR.depth=depth_rw
UC.APR.elev=He
UC.APR.incl=Hd
UC.calcAr_prod_rate_whole_rock()

#plot
fig,ax = plt.subplots()
ax.semilogx(UC.APR.P39Ar_ev_n,depth_r,'r--',label='n$_{ev}$')
ax.semilogx(UC.APR.P39Ar_mu_n,depth_r,'b--',label='n$_{\mu}$')
ax.semilogx(UC.APR.P39Ar_alpha_n,depth_r,'g--',label='n$_{UTh}$')
ax.semilogx(UC.APR.P39Ar_mu,depth_r,'c--',label='$\mu$')
ax.semilogx(UC.APR.Ar39_value,depth_r,'k-',label='Total')
ax.set_ylim(50,0)
ax.set_xlim(10**-4,10**2)
ax.set_xlabel('atoms/g_rock/yr')
ax.set_ylabel('Depth (m)')
ax.legend()
plt.show()
