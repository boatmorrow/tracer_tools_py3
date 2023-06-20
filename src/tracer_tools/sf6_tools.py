#set of utilities for working with dissolved sf6.  
#W. Payton Gardner - University of Montana - Dept. of Geosciences
import tracer_tools.noble_gas_tools as ng
import numpy as N
from tracer_tools.cfc_tools import get_gas_conc

#removing freesteam for now
#def vapor_pressure_atm(T):
#    """P_vapor = vapor_pressure_atm(T) returns the vapor pressure in GPa for pure water
#    at temperature T in C. Uses IAPWS-IF97 and IAPWS-95 formulation using the freesteam
#    libraries http://freesteam.sourceforge.net/"""
#   S = steam_Tx(T+273.15,0);
#    P_vapor = S.p/1.01325e5;    
#   return P_vapor

def vapor_pressure_atm(T):
    """returns the vapor pressure using nobe gas tools module - uses Antione equation.  See
    noble gas tools for more documentation"""
    P_vapor = ng.vapor_pressure(T);  
    P_vapor = P_vapor/0.000101325;   
    return P_vapor
    
def lapse_rate(Elev):
    """P = lapse_rate(Elev).  Rreturns the atmospheric pressure for elevation Elev (m)
    lapse rate taken from Kip Solomon University of Utah. Pressure in GPa.  
    Assumes a pressure of 1 atm at sea level."""
    P = ((1-.0065*Elev/288.15)**5.2561); #Atm  
    return P


def solubility_sf6(T,S=0.0):
    '''Ksf6 = solubility_sf6(T,S,P). Returns the solubility coefficient 
    in mol atm^-1 kg^-1 for sf6 where C(mol/kg) = Ksf6*z_i*(P-Pw) from Bullister 2002.
    Where T is temp in celcius, and S is salinity in parts per thousand'''
    T_k = T+273.15;
    a_1 = -98.7264000;
    a_2 = 142.803;
    a_3 = 38.8746;
    b_1 = 0.0268696;
    b_2 = -0.0334407;
    b_3 = 0.0070843;
    K = N.exp(a_1+a_2*(100/T_k)+a_3*N.log(T_k/100)+S*(b_1+b_2*(T_k/100)+b_3*(T_k/100)**2));
    return K;
    
def equil_conc_sf6(T,z_i,S=0.0,P=1.):
    '''C = equil_conc_sf6(T,z_i,S,P) returns the equilibrium concentration of SF6
    in fmol/kg, for the given atmospheric concentration, for the given temperature (C), salinity (parts per thousand) and atmospheric pressure (atm).  Where C = Ksf*z_i*(P-Pw). z_i is the air mixing
    ratio in parts per trillion volume.  Historical air mixing ratios will need to
    added into this tool kit. See Bullister 2002'''
    P_vapor = vapor_pressure_atm(T);
    Ksf6 = solubility_sf6(T,S);  #mol atm^-1 l^-1
    C = Ksf6*z_i*(P-P_vapor)/.001; #fmol/kg (because z_i in pptv)
    return C

def ce_exc_conc_sf6(Ae,z_i,F=0.0,T=20.,S=0.0,P=1.):
    '''C = ce_exc_conc_sf6(Ae,z_i,T=20.,F=0.0,P=1.) returns the concentration of SF6 
    in fmol/kg due to excess air, for the given temperature (C), salinity (parts per thousand) and 
    atmospheric pressure (atm).  Where C = Ksf*z_i*(P-Pw). z_i is the air mixing 
    ratio in parts per trillion volume.'''
    P_vapor = vapor_pressure_atm(T);
    Ksf6 = solubility_sf6(T,S);  #mol atm^-1 l^-1
    C_eq = Ksf6*z_i*(P-P_vapor)/.001; #fmol/kg (because z_i in pptv)
    Ae = Ae/22414/1e-15*1000; #Ae in
    C_exc = ((1-F)*Ae*(z_i*1e-12))/(1+F*Ae*((z_i*1e-12/C_eq)));
    return C_exc

def equil_air_conc(C_meas,T_rech,E_rech,Ae=0.0,F=0.0,S_rech=0.0,):
    '''z_i = equil_air_conc(C_meas,T_rech,Ae,F) returns the equilibrium air concentration
    given a measured concentation and the estimated recharge temperature in celcius(T_rech), 
    excess air in ccSTP/kg (Ae),  and fractionation factor F (from Aeshbach Hertig 200).  In
    the case of F=0 reverts to the complete dissolution model.'''
    Ki = solubility_sf6(T_rech,S_rech);
    P_da = lapse_rate(E_rech)-vapor_pressure_atm(T_rech);
    molar_volume = 22414.1;
    z_i = (C_meas + ((C_meas*F*(Ae/molar_volume))/(Ki*P_da)))/(Ki*P_da+(Ae/molar_volume))*.001;
    return z_i
    
def age_date(z_i,hemisphere='NH'):
    '''returns the year where the measured air mixing ratio is closest to the value z_i. 
    for the hemisphere of interest'''
    df = get_gas_conc()
    Rech_year = N.idxmin(abs(df['SF6'+hemisphere]-z_i)).year
    return Rech_year
