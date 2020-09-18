# -*- coding: utf-8 -*-
"""
evaporation_models.py
"""
import numpy as N
import pdb
from pylab import *
import matplotlib as mpl
from noble_gas_tools import *



def read_array(filename, dtype, separator='\t'):
    """ Read a file with an arbitrary number of columns.
        The type of data in each column is arbitrary
        It will be cast to the given dtype at runtime
    """
    #pdb.set_trace();
    cast = N.cast
    data = [[] for dummy in range(len(dtype))]
    for line in open(filename, 'r'):
        fields = line.strip().split(separator)
        for i, number in enumerate(fields):
            #print i;
            data[i].append(number)
    for i in range(len(dtype)):
        data[i] = cast[dtype[i]](data[i])
    return N.rec.array(data, dtype=dtype)

def mix(C_1,C_2,nsteps=1000.):
    """C_mix = mix(C_1,C_2,nsteps=100) returns a vector length nsteps+1 starting with a concentration of 
    C_1 and ending with a concentration of C_2 according to mass balance - C_t = xC_1 + (1-x)C_2 where
    x is the mass fraction of water of concentration C_1"""
    step = (1.-0)/nsteps;
    x = arange(1,step,-step);
    C_mix = x*C_1 + (1-x)*C_2;
    return C_mix

##need to set mydescr to read in the latest noble gas data as txt extracted from Yellowstone_data_compilation_8.xls
#mydescr = N.dtype([('Sample Name', 'a5'), ('Temp', 'f4'),('Cl','f4'),('d18O','f4'),('dD','f4')]);
#myrecarray = read_array('stable_isotope_plot_data.txt', mydescr);
#
##extract the cool, warm, and hydrothermal temperature bins
#cool_Cl = myrecarray['Cl'][myrecarray['Temp']<20.];
#cool_18O = myrecarray['d18O'][myrecarray['Temp']<20.];
#cool_D = myrecarray['dD'][myrecarray['Temp']<20.];
#
#warm_Cl = myrecarray['Cl'][(myrecarray['Temp']>=20.) & (myrecarray['Temp']<=50.)];
#warm_18O = myrecarray['d18O'][(myrecarray['Temp']>=20.) & (myrecarray['Temp']<=50.)];
#warm_D = myrecarray['dD'][(myrecarray['Temp']>=20.) & (myrecarray['Temp']<=50.)];
#
#therm_Cl = myrecarray['Cl'][myrecarray['Temp']>50.];
#therm_18O = myrecarray['d18O'][myrecarray['Temp']>50.];
#therm_D = myrecarray['dD'][myrecarray['Temp']>50.];

#mpl.rcParams['text.usetex'] = True;
#mpl.rcParams['font.family'] = 'serif';
#mpl.rcParams['font.serif'] = 'Times, Palatino, New Century Schoolbook, Bookman, Computer Modern Roman';
#mpl.rcParams['font.sans-serif'] = 'Helvetica, Avant Garde, Computer Modern Sans serif';
#mpl.rcParams['font.cursive'] = 'Zapf Chancery';
#mpl.rcParams['font.monospace'] = 'Courier, Computer Modern Typewriter';
mpl.rcParams['font.size'] = 18;
mpl.rcParams['axes.labelsize'] = 18;
mpl.rcParams['xtick.labelsize'] = 16;
mpl.rcParams['ytick.labelsize'] = 16;
mpl.rcParams['axes.titlesize']= 20;
#mpl.rcParams['text.dvipnghack'] = 'False';
mpl.rcParams['figure.figsize'] = 12, 8;


#add in evaporation models
R_h = .9
T = 5.
#pdb.set_trace();
f = N.arange(1,4e-4,-.0001);
delta_o = -15.0
R_o = delta_o/1000+1;
#pdb.set_trace()
alpha_mean = alpha_oxygen(T);
alpha_evap = alpha_mean;
R_a = R_o/(alpha_mean);
delta_eps = 14.2*(1-R_h)/1000;
eps = (alpha_mean - 1);
A = (R_h*R_a + delta_eps + (eps/alpha_evap))/(1 - R_h + delta_eps);
B = (R_h - delta_eps - (eps/alpha_evap))/(1 - R_h + delta_eps);
R_evap = N.ones(len(f));
delta_evap = N.ones(len(f));
for i in range(len(f)):
    R_evap[i] = (R_o - (A/B))*f[i]**B + A/B;
    delta_evap[i] = (R_evap[i]-1)*1000;
Cl_o = .1
Cl_evap = Cl_o/f;
subplot(2,1,1);
semilogx(Cl_evap,delta_evap);
show();
