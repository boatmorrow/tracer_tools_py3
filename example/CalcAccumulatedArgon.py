import tracer_tools.noble_gas_tools as ng
import tracer_tools.He_tools as He
import numpy as np

#create He tools rock type
UC = He.rock_type()

#set to average upper crust as determined from Ballentine and Burnanrd
UC.defined_rock_type() #upper crust is the default and only defined rock type for now

#calcuate Helium production rate
He = 100 #elevation
Hd = 65 #magnetic inclination
depth = 0.5#mwe
depth = 0.5*100 #1000k/m3 * 1/10 kg/m2 -> g/cm2
UC.APR.depth=depth
UC.APR.incl=Hd
UC.APR.elev=He

print('He4 production rate is %1.3g '%UC.HPR.four_He_value, UC.HPR.units)

#calculate Argon production rate
UC.calcAr_prod_rate_whole_rock() #calculates 40Ar production

print('Ar40 production rate is %1.3g '%UC.APR.Ar40_value, UC.APR.units)

#calculate accumulated Argon
age = 300e6  # 300 million years ~ creteaceous granite.
print('Age is %1.1f Myr' %(age/1.e6))
UC.calcAr40_accum(age) #age needs to be in years

UC.switch_units('ccSTP/g_rock/yr')

print('Accumulated Argon is %1.3g'%UC.APR.Ar40_accum, '/'.join(UC.APR.units.split('/')[0:2]))

#desired amount of Argon
V_d = 3e-3 #ccSTP - 3 microliters of Argon

#release fraction V_rel/V_tot
xi = .8 #temperature dependent
C_r = UC.APR.Ar40_accum

#required grams of rock
g_r = V_d/C_r/xi
print('Required grams of rock %1.3f'%g_r, UC.APR.units.split('/')[1])

#required volume of rock 
n = 0.4 #porosity of unconsolidated sediment (conservative)
rho_s = 2.7 # g/cc density of rock grains - alumino-silcate
rho_b = rho_s*(1-n) #g/cc_rev
V_por = g_r/rho_b #cc pulverized rock

print('Required cc of rock %1.3g'%V_por)


#required length for 1cm radius tube
r = 1 #cm radius
A_cyl = np.pi*r**2
h_req = V_por/A_cyl

print('Required length for a tube of radus %1.2f cm is  %1.3f cm' %(r,h_req))

