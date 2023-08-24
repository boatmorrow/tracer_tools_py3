# 07.08.2023

import numpy as np
from tracer_tools.noble_gas_tools import lapse_rate
import pdb
import matplotlib.pyplot as plt
import tracer_tools.He_tools as he

#rock type
VG = he.rock_type()
VG.defined_rock_type(rname='VariscanGranite')

#heidelberg depth profile
He = 0 #elevation
Hd = 80 #magnetic inclination
depth_r = np.linspace(0,50) #m
depth_rw = depth_r*250 #assuming 2500kg/m3 rock density
VG.APR.depth=depth_rw
VG.APR.elev=He
VG.APR.incl=Hd
VG.APR.Ca42=0
VG.calcAr_prod_rate_whole_rock()
VG.calc39ArSecEqul()
print('Equilibrium 39Ar concentration a surface is %1.3g' %VG.APR.Cseq[0])

# Calculating the Age Rock sample Heidelberg Granite, measured in ATTA Juli 10th 2023

#Extraction information:
V = 3.5e-3 #ccSTP extracted total gas volume
frac_Ar = 0.848 #fraction of Argon in the sample
V_Ar = frac_Ar*V  
# uncertainty for Argon volume: ~ 0.5% [from Master's Thesis Philip Hopkins 2018], negligible compared to uncertainty of measured Ar39/Ar ratio 
dilution = 1000/(V_Ar*1e3) #dilution factor (1000muL to 3.5 muL)
m_rock = 220 #[g] based on precision scale +-0.1g, error negligible
X_ex = 0.8 #fraction of argon released from rock, ad hoc assumption
d_X_ex = 0.2 

#ATTA data
ratio_atta = 2.52 #R/R_a
err_low = 0.24
err_high = 0.26 
err_avg = (err_low+err_high)/2 #for now, error propagate assuming symmetric gaussian
#corrected for dilution
R_corr = ratio_atta*dilution #R/R_a undiluted
d_R_corr = err_avg*dilution
R_a = 8e-16 #Benetti et al. 2007
# err_R_a = 0.6e-16 neglected, since R_a is used as a normalisation constant

#constants
ccSTP_to_mol = 1/22414 #mol/ccSTP
N_Av = 6.02214076e23 #mol^-1
molM_Ar40 = 39.962 #g/mol
molM_Ar39 = 38.964 #g/mol
molM_Ka =  39.0983 #g/mol standard atomic weight of potassium
lam39 = np.log(2)/268 #[yr^-1] decay constant of argon, half-life 268+-8 years
d_lam39 = lam39*8/268 #uncertainty


#assumptions: 
r_age = 330e6 #age of the granite
d_r_age = 100e6


#calculating sample Argon 39 concentration in atoms/g_rock based on measured Argon 40 content in the gas extraction 
Nt_extr40 = R_corr * R_a * N_Av * ccSTP_to_mol * V_Ar / X_ex / m_rock #0.8 is the fraction of Argon released
d_Nt_extr40 = Nt_extr40 * np.sqrt((d_R_corr/R_corr)**2 + (d_X_ex/X_ex)**2)
print('39Ar content based on measured 40Ar from gas extraction: %2.2f +- %2.2f atoms/g_rock' %(Nt_extr40, d_Nt_extr40))

#calculation sample Argon 39 concentration based on theoretical Argon 40 content from rock age estimate
# Argon 40 content 
lamK40 =  5.543-10 #[yr^−1] decay constant of K40
lamK40_to_Ar40 = 5.808e-11 #[yr^−1]
K40 = 5.92e-6 * N_Av/molM_Ka #K40 atoms/g_rock, 5.92e.6=weight fraction of K*isotopic ratio of K40/K
d_K40 = 0.08e-6 *N_Av/molM_Ka #based on 1% uncertainty of potassium content and 1/117 in isotope ratio 
VG.calcAr40_accum(r_age) #calculation N40 content with He_tools
N40 = VG.APR.Ar40_accum 
#error of argon 40 content
dN40dK40 = lamK40_to_Ar40/lamK40*(np.exp(lamK40*r_age)-1) #partial derivative with respect to K40 content
dN40dr_age = K40*lamK40_to_Ar40*np.exp(lamK40*r_age) #partial derivative with respect to rock age
d_N40 = np.sqrt(dN40dK40**2 * d_K40**2 + dN40dr_age**2 * d_r_age**2) #adding errors in quadrature

# Argon 39 concentration
Nt_aged40 = R_corr*R_a*N40
d_Nt_aged40 = Nt_aged40*np.sqrt((d_R_corr/R_corr)**2 +(d_N40/N40)**2)
print('39Ar content based on calculated 40Ar: %2.2f +- %2.2f atoms/g_rock' %(Nt_aged40, d_Nt_aged40))

# initial concentration
N0 = VG.APR.P39Ar_alpha_n[-1]/lam39
d_N0 = 0.3*N0 # estimate taken from Sramek 2017, uncertainty of lamba is negligible
print('initial concentration at first exposure is assumed to be %2.2f +- %2.2f atoms/g_rock' %(N0, d_N0))
# secular equilibrium
Ceq = VG.APR.Cseq[0]
d_Ceq = 0.5*Ceq # estimate based on Musy 2023
print('Surface secular equilibrium 39Ar concentration is %2.2f +- %2.2f atoms/g_rock' %(Ceq, d_Ceq))

#exposure age
def exp_age(lamma, Nt, eqCon, N_i):
    age = 1/lamma * np.log((N_i-eqCon)/(Nt-eqCon))
    return age
def err_age(lamma, err_lam, Nt, err_Nt, eqCon, err_eqCon, Ni, err_Ni): #for now not including errors of the production rate
    dtdlam = exp_age(lamma, Nt,eqCon, Ni)/lamma
    dtdN_i = 1/(Ni*lamma - eqCon*lamma)
    dtdNt = 1/(eqCon*lamma - Nt*lamma)
    dtdeqCon = (Nt-Ni)/lamma/(Ni-eqCon)/(eqCon-Nt)
    err_2 = dtdlam**2 * err_lam**2 + dtdN_i**2 * err_Ni**2 + dtdNt**2 * err_Nt**2 + dtdeqCon**2 * err_eqCon**2
    err = np.sqrt(err_2)
    return err

age_extr40 = exp_age(lam39, Nt_extr40, Ceq, N0)
d_age_extr40 = err_age(lam39, d_lam39, Nt_extr40, d_Nt_extr40, Ceq, d_Ceq, N0, d_N0)
print('Rock exposure age based on measured 40Ar content: %4.2f +- %3.2f yrs' %(age_extr40, d_age_extr40) )

age_40theo = exp_age(lam39, Nt_aged40, Ceq, N0)
err_age_40theo = err_age(lam39, d_lam39, Nt_aged40,  d_Nt_aged40, Ceq, d_Ceq, N0, d_N0)
print('Rock exposure age based on calculated 40Ar content: %4.2f +- %3.2f yrs' %(age_40theo, err_age_40theo) )    


#measurement range 
N_array = np.arange(N0, Ceq, 5)
plt.figure(figsize=(10,8))
plt.plot(N_array, exp_age(lam39, N_array, Ceq, N0))
plt.xlabel('Ar39 content [atoms/g_rock]')
plt.ylabel('Surface exposure age [yr]')
#plt.yscale('log')
plt.grid(True)

# uncertainty: dt/dN = [(Ceq-N)*lamma]^-1
N_array = np.arange(N0+3, 4500, 5)
t = exp_age(lam39, N_array, Ceq, N0)
def del_age(delN_rel, N, eqCon, lamma):
    dt = delN_rel*N/((eqCon-N)*lamma)
    return dt
plt.figure(figsize=(10,6))
delN = np.array([0.1,0.2,0.5,0.7,1])
for i in range(len(delN)):
    plt.plot(N_array, del_age(delN[i], N_array, Ceq, lam39)/t, label='relative delta_N %2.1f' %(delN[i]))
plt.grid(True)
plt.legend()
plt.xlabel('Ar39 concentration N [atoms/g_rock]')
plt.ylabel('Relative Age uncertainty Delta t/t =dt/dN * Delta N / t')
plt.yscale('log')
plt.show()