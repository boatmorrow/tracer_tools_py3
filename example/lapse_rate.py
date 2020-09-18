# -*- coding: utf-8 -*-
"""
Cacluate noble gas concentrations for an elevation increase and plot them up.  Assumes an average pressure and temperature lapse rate.
"""

from noble_gas_tools import *
import numpy as N
from pylab import *

t_lr = 5. #5C/km a typical wet lapserate
t_lr_dry = 9.8 #the dry adiabatic lapse rate

E = N.arange(0,1000)
T = 10 - t_lr/1000.*E
P = []
Ne = []
He = []
Xe = []
Ar = []
for i in range(len(E)):
    Pi = lapse_rate(E[i])
    #Pi = lapse_rate(0)
    P.append(Pi)
    Ne.append(equil_conc('Ne',T[i],P=Pi))
    He.append(equil_conc('He',T[i],P=Pi))
    Xe.append(equil_conc('Xe',T[i],P=Pi))
    Ar.append(equil_conc('Ar',T[i],P=Pi))


figure()

subplot(2,1,1)
plot(E,T,label='temp C')
ylabel('T (C)')

subplot(2,1,2)
plot(E,P,label='Pres. GPa')
ylabel('P (GPa)')

figure()
subplot(2,2,2)
plot(E,Ne,label='Ne ccSTP/g')
#xlabel('Elev. (m)')
ylabel('[Ne]$_aq$ (ccSTP/g)')

subplot(2,2,4)
plot(E,Xe,label='Xe ccSTP/g')
#xlabel('Elev. (m)')
ylabel('[Xe]$_aq$ (ccSTP/g)')

subplot(2,2,1)
plot(E,He,label='He ccSTP/g')
xlabel('Elev. (m)')
ylabel('[He]$_aq$ (ccSTP/g)')

subplot(2,2,3)
plot(E,Ar,label='Ar ccSTP/g')
xlabel('Elev. (m)')
ylabel('[Ar]$_aq$ (ccSTP/g)')
show()
