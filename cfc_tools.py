#suite of tools for working with disolved cfcs in groundwater
#W. Payton Gardner - University of Montana, Dept. of Geosciences.

#from freesteam import *
import noble_gas_tools as ng
import numpy as N
import datetime
import pandas as pd
import pdb
import subprocess


def update_gas_conc(url='http://cdiac.ess-dive.lbl.gov/ftp/oceans/CFC_ATM_Hist/CFC_ATM_Hist_2015/CFC_atmospheric_histories_revised_2015_Table1.csv'):
    '''Grab the latest and greates sf6 and cfc gas concentrations.  The URL will need to be updated periodically.  Data from:
        Bullister, J.L. 2015. Atmospheric Histories (1765-2015) for CFC-11, CFC-12, CFC-113, CCl4, SF6 and N2O. NDP-095(2015). http://cdiac.ornl.gov/ftp/oceans/CFC_ATM_Hist/CFC_ATM_Hist_2015. Carbon Dioxide Information Analysis Center, Oak Ridge National Laboratory, US Department of Energy, Oak Ridge, Tennessee. doi: 10.3334/CDIAC/otg.CFC_ATM_Hist_2015.
        '''
    subprocess.call(['curl','-o','cfc_sf6_atm.txt',url])

def get_gas_conc(f='cfc_sf6_atm.txt'):
    df = pd.read_csv(f,header=0)
    units = df.loc[0]
    df = df.drop(0)
    #Fractional Years!
    ix = pd.date_range(datetime.datetime(int(df.loc[df.index[0]]['Year']),1,1),datetime.datetime(int(df.loc[df.index[-1]]['Year']),1,1),freq='AS')
    dt = pd.Timedelta(6,'M')
    ix = ix+dt
    df = df.set_index(ix)
    df = df.astype(N.float64)
    return df

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

def solubility_cfc(T,cfc=12,S=0.0):
    '''Kcfc11 = solubility_cfc11(T,S). Returns the solubility coefficient 
    in mol atm^-1 kg^-1 for cfc where C(mol/kg) = Ksf6*z_i*(P-Pw) from Bullister 2002.
    Where T is temp in celcius, and S is salinity in parts per thousand'''
    T_k = T+273.15;
    #pdb.set_trace();
    cfc_dict = {11:[-136.2685,206.1150,57.2805,-.148598,0.095114,-0.0163396],12:[-124.4395,185.4299,51.6383,-0.149779,0.094668,-0.0160043],113:[-136.129,206.475,55.8957,-0.02754,0.006033,0,]};
    a_1 = cfc_dict[cfc][0];
    a_2 = cfc_dict[cfc][1];
    a_3 = cfc_dict[cfc][2];
    b_1 = cfc_dict[cfc][3];
    b_2 = cfc_dict[cfc][4];
    b_3 = cfc_dict[cfc][5];
    K = N.exp(a_1+a_2*(100/T_k)+a_3*N.log(T_k/100)+S*(b_1+b_2*(T_k/100)+b_3*(T_k/100)**2));
    return K;
    
def equil_conc_cfc(T,z_i,cfc=12,S=0.0,P=1.):
    '''C = equil_conc_cfc(T,z,cfc=12,S=0,P=1.) returns the equilibrium concentration of CFC 
    in pmol/kg, for the given temperature (C), zir mixing ration (pptv), salinity (parts per thousand) and 
    atmospheric pressure (atm).  Where C = Ksf*z_i*(P-Pw). z_i is the air mixing 
    ratio in parts per trillion.  Historical air mixing ratios will need to
    added into this tool kit. See Bullister 2002'''
    P_vapor = vapor_pressure_atm(T);
    Kcfc = solubility_cfc(T,cfc,S);  #mol atm^-1 l^-1
    C = Kcfc*z_i*(P-P_vapor); #pmol/kg (because z_i in pptv)
    return C

def equil_air_conc(C_meas,T_rech,E_rech,cfc=12,Ae=0.0,F=0.0,S_rech=0.0):
    '''z_i = equil_air_conc(C_meas,T_rech,Ae,F) returns the equilibrium air concentration
    given a measured concentation and the estimated recharge temperature in celcius(T_rech), 
    excess air in ccSTP/kg (Ae),  and fractionation factor F (from Aeshbach Hertig 200).  In
    the case of F=0 reverts to the complete dissolution model.'''
    Ki = solubility_cfc(T_rech,cfc,S_rech);
    P_da = lapse_rate(E_rech)-vapor_pressure_atm(T_rech);
    molar_volume = 22414.1;
    z_i = (C_meas + ((C_meas*F*(Ae/molar_volume))/(Ki*P_da)))/(Ki*P_da+(Ae/molar_volume));
    return z_i
    
def age_date(z_i,cfc=12,hemisphere='NH'):
    '''returns the year where the measured air mixing ratio is closest to the value z_i. 
    Right now for northern hemisphere only, but adding southern hemisphere is easy peasy'''
    name_dict = {11:'CFC11'+hemisphere,12:'CFC12'+hemisphere,113:'CFC113'+hemisphere}
    #dd = N.load('cfc_data.dat.npy');
    df = get_gas_conc()
    Rech_year = N.argmin(abs(df[name_dict[cfc]]-z_i)).year
    return Rech_year
    
def rayleigh_frac_cfc(T_sep,T_ref,S=0):
    """ CFC11_res, CFC12_res = rayleight_frac_cfc(T_sep,T_ref,S=0). Returns residual concentrations 
    of CFC11 and 12 for a rayleigh type distillation.  Where T_sep is the temperature of 
    separation, and T_ref is the initial temperature of atmospheric equilibration, and S is the salinity
    at separation (assumes zero salinity at recharge) in parts per thousand.  Right now initial concentrations
    of CFC are assumed to be from year 1991.  Equation taken from 
    Porcelli, D.; Ballentine, C. J. & Wieler, R. (ed.) Noble Gases in Geochemistry and 
    Cosmochemistry Mineralogical Society of America , 2002, 47 pg 554. Where T_ref is the temperature
    of ASW that the sample is referenced too and T_sep is the temperature of vapor phase separation in 
    celcius,S is salinity in mol/l?, and CFC11_frac is the residual cfc11."""
    CFC11_frac = N.arange(1,.001,-.001);
    #pdb.set_trace();
    K_i = 1/solubility_cfc(T_sep,cfc=12);
    K_Ar = 1/solubility_cfc(T_sep,cfc=11);
    alpha = K_i/K_Ar;
    i_0 = equil_conc_cfc(T_ref,496,cfc=12,S=S);
    Ar_0 = equil_conc_cfc(T_ref,270.4,cfc=11,S=S);
    frac_0 = i_0/Ar_0;
    frac_i = frac_0*CFC11_frac**(alpha-1);
    CFC11_res = CFC11_frac*Ar_0;
    CFC12_res = CFC11_res * frac_i;
    return CFC11_res, CFC12_res
