import tracer_tools.noble_gas_tools as ng
import tracer_tools.He_tools as He
import numpy as np
import matplotlib.pyplot as mpl





#inputs
#desired amount of Argon
#V_d = 1e-4 #ccSTP - 0.1 microliters of Argon
V_d = 3e-3 #ccSTP - 3 microliters of Argon

#release fraction V_rel/V_tot
xi = .7 #temperature dependent

#required volume of rock 
n = 0.4 #porosity of unconsolidated sediment (conservative)
rho_s = 2.7 # g/cc density of rock grains - alumino-silcate

#extraction finger radius
r = 1 #cm radius



#create He tools rock type
UC = He.rock_type()

#set to average upper crust as determined from Ballentine and Burnanrd
UC.defined_rock_type() #upper crust is the default and only defined rock type for now

#calcuate Helium production rate
UC.calcHe_prod_rate_whole_rock() #calculates 4He and 3He production

print('He4 production rate is %1.3g '%UC.HPR.four_He_value, UC.HPR.units)

#calculate Argon production rate
UC.n_sf_spec = 'C:/Users/Fanka Neumann/Documents/Codes_MA/tracer_tools_py3/data/JamtalNeutronFlux.CSV'
UC.calcAr_prod_rate_whole_rock() #calculates 40Ar production

print('Ar40 production rate is %1.3g '%UC.APR.Ar40_value, UC.APR.units)

#calculate accumulated Argon
age = np.logspace(7,9)  # 300 million years ~ creteaceous granite.
UC.calcAr40_accum(age) #age needs to be in years
UC.switch_units('ccSTP/g_rock/yr')
C_r = UC.APR.Ar40_accum


#required grams of rock
g_r = V_d/C_r/xi


#porous media calcs
rho_b = rho_s*(1-n) #g/cc_rev
V_por = g_r/rho_b #cc pulverized rock

#required length for tube
A_cyl = np.pi*r**2
h_req = V_por/A_cyl


#make plots 
fig1,ax1 = mpl.subplots(figsize=(10,8))
ax1.loglog(age/1e6,C_r,'r-',lw=2)
ax1.set_xlabel('Rock Age (Mya)')
ax1.set_ylabel('C$^{40}$Ar$_{rock}$  ' + '/'.join(UC.APR.units.split('/')[0:2]))


fig2,axes = mpl.subplots(nrows=2,figsize=(10,8))
axes[0].loglog(age/1e6,V_por,'k-',lw=2)
axes[0].set_xlabel('Rock Age (Mya)')
axes[0].set_ylabel('Volume Crushed Rock (cc)')
axes[0].vlines(250,10,1000,color='k',alpha=0.7,lw=2,label='Average Crust')
axes[0].grid(True)


#axes[1].loglog(age/1e6,h_req,'k-',lw=2)
#axes[1].set_xlabel('Rock Age (Mya)')
#axes[1].set_ylabel('Height of rock (r=%1.1f) (cm)'%r)

axes[1].plot(age/1e6,g_r,'k-',lw=2)
axes[1].set_xlabel('Rock Age (Mya)')
axes[1].set_xlim(240,260)
axes[1].set_ylabel('Required mass (g)')
axes[1].set_ylim(100,300)
axes[1].vlines(250,10,1000,color='k',alpha=0.7,lw=2,label='Average Crust')
axes[1].grid(True)
#mpl.savefig("/Users/wpgardner/pres/2023_Goldschmidt/RequiredMass.png",dpi=200)
#ax1.fill_between()
mpl.show()

#CALCULATE      