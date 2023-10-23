import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.interpolate import interp1d

 #read in the cross section
inp_file = 'C:/Users/Fanka Neumann/Documents/Codes_MA/tracer_tools_py3/src/tracer_tools/K39_XS.csv'
K39xsec = pd.read_csv(inp_file)
inp_file = 'C:/Users/Fanka Neumann/Documents/Codes_MA/tracer_tools_py3/src/tracer_tools/Ca42_XS.csv'
Ca42xsec = pd.read_csv(inp_file)

 #read in cosmic spectral flux from excel sheet for heidelberg lattitude and elevation
inp_file = "C:/Users/Fanka Neumann/Documents/Codes_MA/tracer_tools_py3/data/HeidelbergNeutronFlux.csv"
dfsc = pd.read_csv(inp_file) #neutrons/cm^2/s/MeV
#truncate to the cross section spectrum
dfstc = dfsc[dfsc['Energy']<K39xsec['Energy [MeV]'].max()] #for K
#fit a function to xsec
f = interp1d(K39xsec['Energy [MeV]'],K39xsec['Cross-Section [barns]']*1e-24)
#interpolate to neutron flux
sigma_newc=f(dfstc['Energy'])
#evaporation
intgrnd_ev=sigma_newc*dfstc['Neutron'] #for a particular depth



print(K39xsec.head())

plt.figure(figsize=(10,6))
plt.plot(K39xsec['Energy [MeV]'], K39xsec['Cross-Section [barns]'], label='K39(n,p)Ar39')
plt.plot(Ca42xsec['Energy [MeV]'], Ca42xsec['Cross-Section [barns]'], label='Ca42(n,alpha)Ar39')
plt.xlabel('Energy [MeV]')
plt.ylabel('Reaction cross-section [barns]')
#plt.yscale('log')
plt.ylim(0,0.4)
#plt.xscale('log')
plt.xlim(0,30)
plt.grid(True)
plt.legend()
plt.title('Reaction cross-sections for Argon-39 production')
# plt.savefig('C:/Users/Fanka Neumann/Documents/Codes_MA/Cross_sections_linear_plot.pdf')
plt.show()

plt.figure(figsize=(10,6))
plt.plot(dfstc['Energy'], intgrnd_ev)
plt.xlabel('Energy [MeV]')
plt.xscale('log')
plt.ylabel('Cross-section * surface particle flux [neutrons/s/MeV]')
plt.yscale('log')
plt.title('Integrand for Argon production')
plt.grid(True)
plt.show()