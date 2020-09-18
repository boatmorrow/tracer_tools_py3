# -*- coding: utf-8 -*-
"""
Created on Wed Jul 06 08:38:09 2011

@author: gar279
"""

from noble_gas_tools import *
from numpy import *
from pylab import *

T = arange(10,40);
ec_v = vectorize(equil_conc);
C_T = ec_v('Ne',T);
C_T_Xe = ec_v('Xe',T);
C_T_Ar = ec_v('Ar',T);
C_T_Kr = ec_v('Kr',T);
plot(T,C_T/min(C_T));
plot(T,C_T_Ar/min(C_T_Ar));
plot(T,C_T_Kr/min(C_T_Kr));
plot(T,C_T_Xe/min(C_T_Xe));
ylabel('C/C$_{min}$');
xlabel('T $^\circ$C');
legend(('Neon','Argon','Krypton','Xenon'));

figure();
F=0;
Ae_v = linspace(0,.1);
ce_exc_v = vectorize(ce_exc);
C_ae = ce_exc_v('Ne',F,Ae_v,25.) + equil_conc('Ne',25);
C_ae_ar = ce_exc_v('Ar',F,Ae_v,25.) + equil_conc('Ar',25);
C_ae_kr = ce_exc_v('Kr',F,Ae_v,25.) + equil_conc('Kr',25);
C_ae_xe = ce_exc_v('Xe',F,Ae_v,25.) + equil_conc('Xe',25);
plot(Ae_v,C_ae/min(C_ae));
plot(Ae_v,C_ae_ar/min(C_ae_ar));
plot(Ae_v,C_ae_kr/min(C_ae_kr));
plot(Ae_v,C_ae_xe/min(C_ae_xe));
xlabel('Ae (ccSTP/g)');
ylabel('C/C$_{min}$');
legend(('Neon','Argon','Krypton','Xenon'));


figure();
Ae = .02
F_v=linspace(0,1.);
C_ae_F = ce_exc_v('Ne',F_v,Ae,25.) + equil_conc('Ne',25);
C_ae_ar_F = ce_exc_v('Ar',F_v,Ae,25.) + equil_conc('Ar',25);
C_ae_kr_F = ce_exc_v('Kr',F_v,Ae,25.) + equil_conc('Kr',25);
C_ae_xe_F = ce_exc_v('Xe',F_v,Ae,25.) + equil_conc('Xe',25);
plot(F_v,C_ae_F/min(C_ae_F));
plot(F_v,C_ae_ar_F/min(C_ae_ar_F));
plot(F_v,C_ae_kr_F/min(C_ae_kr_F));
plot(F_v,C_ae_xe_F/min(C_ae_xe_F));
xlabel('F');
ylabel('C/C$_{min}$');
legend(('Neon','Argon','Krypton','Xenon'));
show();
