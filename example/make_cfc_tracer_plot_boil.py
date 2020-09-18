# -*- coding: utf-8 -*-
"""
Created on Tue Jan 04 08:03:56 2011

@author: gar279
"""

import numpy as N
from pylab import *
from cfc_tools import *

T_atm = 4.;
P_atm = lapse_rate(2500);

#load cfc atmospheric concentrations
dd = N.load('cfc_data.dat.npy');

#need to set mydescr to read in the latest noble gas data as txt extracted from Yellowstone_data_compilation_8.xls
#mydescr = N.dtype([('Sample Name', 'a5'), ('Temp', 'f4'), ('trit', 'f4'), ('cfc12', 'f4'),('cfc11','f4')]);
#myrecarray = read_array('tracer_plot_data_cfc11.txt', mydescr);

myrecarray = loadtxt('tracer_plot_data_cfc11.txt', dtype = {'names':('Sample Name','Temp','trit','cfc12','cfc11'),'formats':('a5','f4','f4','f4','f4')});

#extract the cool, warm, and hydrothermal temperature bins
cool_trit = myrecarray['trit'][myrecarray['Temp']<20.];
cool_cfc12 = myrecarray['cfc12'][myrecarray['Temp']<20.];
cool_cfc11 = myrecarray['cfc11'][myrecarray['Temp']<20.];

warm_trit = myrecarray['trit'][(myrecarray['Temp']>=20.) & (myrecarray['Temp']<=50.)];
warm_cfc12 = myrecarray['cfc12'][(myrecarray['Temp']>=20.) & (myrecarray['Temp']<=50.)];
warm_cfc11 = myrecarray['cfc11'][(myrecarray['Temp']>=20.) & (myrecarray['Temp']<=50.)];

therm_trit = myrecarray['trit'][myrecarray['Temp']>50.];
therm_cfc12 = myrecarray['cfc12'][myrecarray['Temp']>50.];
therm_cfc11 = myrecarray['cfc11'][myrecarray['Temp']>50.];


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
cfc11_aq = equil_conc_cfc(T_atm,dd['cfc_11'],cfc=11);
cfc12_aq = equil_conc_cfc(4.,dd['cfc_12']);

# get boiled concentrations
cfc11_res, cfc12_res = rayleigh_frac_cfc(90,4);

#lets have a looksie
figure();
plot(cfc11_aq,cfc12_aq);
plot(cfc11_res,cfc12_res,'--');
plot(cool_cfc11,cool_cfc12,'b^',label='cool');
plot(warm_cfc11,warm_cfc12,'gs',label='warm');
plot(therm_cfc11,therm_cfc12,'ro',label='hydrotherm.');
xlabel(r'[CFC-11]$_{aq}$ (pmol/kg)');
ylabel(r'[CFC-12]$_{aq}$ (pmol/kg)');