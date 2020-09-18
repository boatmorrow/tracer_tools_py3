# -*- coding: utf-8 -*-
"""
Created on Mon Feb 25 08:51:03 2013
make atmospheric input plots
figure 1 will be atmospheric concentration
figure 2 will be concentration in groundwater
@author: wpgardn
"""
from pylab import *
from MakeTotalTimeSeries import *
import matplotlib as mpl


mpl.rcParams['text.usetex'] = False;
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
mpl.rcParams['figure.figsize'] = 14, 10;



ddff = MakeTotalTimeSeries();
ddff_aq = convert2aqueous(ddff,10.0);
ddff_fin = convert2molality(ddff_aq);

figure();
semilogy(ddff.index,ddff['cfc12'],label='CFC-12',linewidth=2);
plot(ddff.index,ddff['cfc11'],label='CFC-11',linewidth=2);
plot(ddff.index,ddff['cfc113'],label='CFC-113',linewidth=2);
plot(ddff.index,ddff['SF6 (ppt)'],label=r'SF$_6$',linewidth=2);
plot(ddff.index,ddff['trit TU'],label=r'$^3$H',linewidth=2);
xlabel('Year');
ylabel('Atmospheric Concentration');
mpl.pyplot.xticks(rotation=45);
legend(loc='best');

figure();
semilogy(ddff_fin.index,ddff_fin['CFC12'],label='CFC-12',linewidth=2);
plot(ddff_fin.index,ddff_fin['CFC11'],label='CFC-11',linewidth=2);
plot(ddff_fin.index,ddff_fin['CFC113'],label='CFC-113',linewidth=2);
plot(ddff_fin.index,ddff_fin['SF6'],label=r'SF$_6$',linewidth=2);
plot(ddff_fin.index,ddff_fin['H3'],label=r'$^3$H',linewidth=2);
xlabel('Year');
ylabel('Recharge Concentration [M]');
mpl.pyplot.xticks(rotation=45);
legend(loc='best');
show();