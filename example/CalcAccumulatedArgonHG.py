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
UC.calcAr_prod_rate_whole_rock() 

print('Ar40 production rate is %1.3g '%UC.APR.Ar40_value, UC.APR.units)

#calculate accumulated Argon
age = 300e6  # 300 million years ~ creteaceous granite.
print('Age is %1.1f Myr' %(age/1.e6))
UC.calcAr40_accum(age) #age needs to be in years

UC.switch_units('ccSTP/g_rock/yr')

print('Theoretical Argon concentration is %1.3g'%UC.APR.Ar40_accum, '/'.join(UC.APR.units.split('/')[0:2]))

#Theoretical amount of Argon in 100g
#C_r = UC.APR.Ar40_accum
#print('Theoretical concentration in rock is %1.3g '%C_r,'/'.join(UC.APR.units.split('/')[0:1]))
g_rock_ext = 100 #grams of rock we extracted from
V_ext_th = UC.APR.Ar40_accum*g_rock_ext
print('Expected volume of 40Ar in 100g is %1.3g muL'%(V_ext_th*1e3))

#release fraction 
V_ext = 2.5 #microliters extracted
xi = V_ext/(V_ext_th*1e3) #temperature dependent
print('Extraction of theoretical fraction after 3 hours is %1.3f'%xi)

V_ext = 3.5 #microliters extracted
xi = V_ext/(V_ext_th*1e3) #temperature dependent
print('Extraction of theoretical fraction after 3 days is %1.3f'%xi)
