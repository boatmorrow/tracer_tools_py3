# -*- coding: utf-8 -*-
"""
Created on Tue Jan 04 08:03:56 2011

@author: gar279
"""

import numpy as N
from pylab import *
from cfc_tools import *
from sf6_tools import *
from He_tools import *
from noble_gas_tools import equil_conc,mix
import pdb

T_atm = 20.;
P_atm = lapse_rate(200)/0.000101325;
#pdb.set_trace();
#load cfc atmospheric concentrations
dd_cfc = N.load('cfc_data.dat.npy');

dd = N.load('sf6_data.dat.npy');


mpl.rcParams['text.usetex'] = True;
#mpl.rcParams['font.family'] = 'serif';
#mpl.rcParams['font.serif'] = 'Times, Palatino, New Century Schoolbook, Bookman, Computer Modern Roman';
#mpl.rcParams['font.sans-serif'] = 'Helvetica, Avant Garde, Computer Modern Sans serif';
#mpl.rcParams['font.cursive'] = 'Zapf Chancery';
#mpl.rcParams['font.monospace'] = 'Courier, Computer Modern Typewriter';
mpl.rcParams['font.size'] = 14;
mpl.rcParams['axes.labelsize'] = 18;
mpl.rcParams['xtick.labelsize'] = 16;
mpl.rcParams['ytick.labelsize'] = 16;
mpl.rcParams['axes.titlesize']= 20;
mpl.rcParams['text.dvipnghack'] = 'False';
mpl.rcParams['figure.figsize'] = 12, 8;



# aqueous concentrations
cfc11_aq = equil_conc_cfc(T_atm,dd_cfc['cfc_11'],cfc=11,P=P_atm);
cfc12_aq = equil_conc_cfc(T_atm,dd_cfc['cfc_12'],P=P_atm);
SF6_aq = equil_conc_sf6(T_atm,dd['conc'],P=P_atm);


#lets have a looksie
figure(1);
plot(cfc12_aq,cfc11_aq,label='Piston Flow');
ylabel(r'[CFC-11]$_{aq}$ (pmol/kg)');
xlabel(r'[CFC-12]$_{aq}$ (pmol/kg)');

# get cfcs for all years in sf6
def find_common_time(dd, cfc12_aq, dd_cfc):
    cfc_ct = N.ones(len(dd['conc']));
    for y in range(len(dd['year'])):
        cfc_ct[y] = cfc12_aq[N.argmin(N.abs(dd['year'][y]-dd_cfc['year']))];
    return cfc_ct


cfc_common_time = find_common_time(dd,cfc12_aq,dd_cfc);


figure(2);
plot(cfc_common_time,SF6_aq);
xlabel(r'[CFC-12]$_{aq}$ (pmol/kg)');
ylabel(r'[SF$_6$]$_{aq}$ (fmol/kg)');

crust_t = define_rock_type(); #average rock type 
crustal_prod = calcHe_prod_rate_whole_rock(crust_t,pth=0.08);
switch_units(crustal_prod,units_desired='ccSTP/g_h2o/yr');
sf6_tau = N.max(dd['year'])-dd['year'];
sf6_tau = sf6_tau[::-1];
tau_vec_long = N.arange(N.ceil(N.max(sf6_tau)),1e6); #years from the most recent sf6 measurment 2006.5
tau_vec = N.concatenate((sf6_tau,tau_vec_long));
sf6_padded = N.concatenate((SF6_aq[::-1],N.zeros(len(tau_vec)-len(dd['year']))));
he4_accum = tau_vec*crustal_prod.four_He_value;
he4_t = he4_accum+equil_conc('He',T_atm,P=P_atm*0.000101325);
he_mix = mix(N.max(he4_accum),equil_conc('He',T_atm,P=P_atm*0.000101325))
sf6_mix = mix(0,N.max(SF6_aq));

figure(3);
semilogy(sf6_padded,he4_t,label='piston flow');
semilogy(sf6_mix,he_mix,label='binary mixture');
xlabel(r'[SF$_6$]$_{aq}$ (fmol/kg)');
ylabel(r'$^4$He$_{aq}$ (ccSTP/g)');
xlim(-.1,1.8)

show();
