# -*- coding: utf-8 -*-
#This is a library of noble gas functions written by Payton Gardner - University of Montana, Dept. of Geosciences
from numpy import exp, log, array, arange, linspace, append
import xlrd
import tkinter
import tkinter.filedialog
import pickle as pickle
#from scipy import interpolate
import pdb
# the atmospheric standard taken from Porcelli et al 2002
atm_std = { 'N2'   : 0.781,
            'O2'   : 0.209,
            'Ar'   : 9.34e-3,
            'CO2'  : 3.7e-4,
            'Ne'   : 1.818e-5,
            'He'   : 5.24e-6,
            'CH4'  : 1.5e-6,
            'Kr'   : 1.14e-6,
            'H2'   : 7e-7,
            'N2O'  : 3e-7,
            'CO'   : 1e-7,
            'Xe'   : 8.17e-8,
            'Rn'   : 6e-20,
            'He3'  : 5.24e-6*0.000140/100,
            'He4'  : 5.24e-6,
            'Ne20' : 1.818e-5*90.50/100,
            'Ne21' : 1.818e-5*0.268/100,
            'Ne22' : 1.818e-5*9.23/100,
            'Ar36' : 9.34e-3*0.3364/100,
            'Ar38' : 9.34e-3*0.0632/100,
            'Ar40' : 9.34e-3*99.60/100,
            'Kr78' : 1.14e-6*0.3469/100,
            'Kr80' : 1.14e-6*2.2571/100,
            'Kr82' : 1.14e-6*11.523/100,
            'Kr83' : 1.14e-6*11.477/100,
            'Kr84' : 1.14e-6*57.00/100,
            'Kr86' : 1.14e-6*17.398/100,
            'Xe124': 8.17e-8*0.0951/100,
            'Xe126': 8.17e-8*0.0887/100,
            'Xe128': 8.17e-8*1.919/100,
            'Xe129': 8.17e-8*26.44/100,
            'Xe130': 8.17e-8*4.070/100,
            'Xe131': 8.17e-8*21.22/100,
            'Xe132': 8.17e-8*26.89/100,
            'Xe134': 8.17e-8*10.430/100,
            'Xe136': 8.17e-8*8.857/100}  

def get_atm_mol_frac(gas):
    '''return the standard mole fraction for the gas and/or isotope of interest.  Ratios taken from Porcelli et al 2002.'''
    pp_i = atm_std[gas]
    return pp_i

def solubility(gas,T,S=0):
    """
    K = solubility(gas,T,S).  Returns the solubility coefficient for the noble gas of interest
    in water with units GPa the satisfies the equation p_i = K_i*x_i where x_i is the mol fraction 
    of gas i in water, K_i is the salt corrected Henry's coefficient and p_i is the partial pressure 
    of gas i in GPa.  If temperature is less the 65 C salting coefficients are included. gas is 
    "He", "Ne", "Ar", "Kr", or "Xe". T is the temperature in degrees C.  S is the salinity in 
    mol l^-1. Solubility and salting coefficients are calculated using date from 
    Porcelli, D.; Ballentine, C. J. & Wieler, R. (ed.) Noble Gases in Geochemistry and 
    Cosmochemistry Mineralogical Society of America , 2002, 47 pg 547
    """
    gas=gas[0:2]
    #pdb.set_trace();
    T_k = T + 273.15;
    # noble gas solubility coefficients
    sol_dict = {'He':[-0.00953,0.107722,0.001969,-0.043825],
                'Ne':[-7.259,6.95,-1.3826,0.0538],
                'Ar':[-9.52,8.83,-1.8959,0.0698],
                'Kr':[-6.292,5.612,-0.8881,-0.0458],
                'Xe':[-3.902,2.439,0.3863,-0.221]};
    #need a Stechnow coefficient dictionary for the other noble gases as well 
    # taken from ng book pg 544 from kennedy
    setch_dict = {'He':[-10.081,15.1068,4.8127],
                'Ne':[-11.9556,18.4062,5.5464],
                'Ar':[-10.6951,16.7513,4.9551],
                'Kr':[-9.9787,15.7619,4.6181],
                'Xe':[-14.5524,22.5255,6.7513]};
    G = setch_dict[gas];
    G1 = G[0];
    G2 = G[1];
    G3 = G[2];
    if T < 65.: #add Setchenow coefficient. note: Setchenow coefficicient are only valid up to 65C
        setch = G1 + (G2/(.01*T_k))+ (G3*log(.01*T_k)); 
        gamma = exp(S*setch);
        #print 'stechenow salinity units uknown need Kennedy 1982 to get a better understanding...'
    else:
        gamma = 1.;
    
    #now calculate the solubility
    A = sol_dict[gas];
    A0 = A[0]  # taken from ng book pg 545 from Crovetto 1981
    A1 = A[1];
    A2= A[2];
    A3 = A[3];
    if gas == 'He':
        ln_F_He = A0 + (A1/(.001*T_k)) + (A2/(.001*T_k)**2) + (A3/(.001*T_k)**3)
        F = exp(ln_F_He);
        Frac_He_gas = 5.24e-6/9.31e-3;
        X_Ar_water = 1./(exp(sol_dict['Ar'][0] + (sol_dict['Ar'][1]/(.001*T_k))+(sol_dict['Ar'][2]/(.001*T_k)**2) + (sol_dict['Ar'][3]/(.001*T_k)**3)))*9.31e-3;
        X_he_water = F*Frac_He_gas*X_Ar_water
        K_h = 5.24e-6/X_he_water;
        K = gamma*K_h;
        
    else:
        ln_K = A0 + (A1/(.001*T_k))+(A2/(.001*T_k)**2) + (A3/(.001*T_k)**3);
        K_h = exp(ln_K); # Henry's coefficient in GPa
        K = gamma*K_h;  # Salinty effects added solubility coeffient in GPa
    return K;

def schmidt(gas,T):
    '''calculate the Schmidt number, the ratio of the kinematic viscosity of water to the diffusion coefficient (for gas exchange)
    for the given gas, at temperature T in celcius.  Schmidt params taken from Raymond 2012.'''
    sch_coeff_dict = {'He': [377.,- 19.15,0.501,- 0.0057],\
                     'Rn': [2939.,- 173.87,4.532,- 0.0468],\
                     'SF6':[3255.,- 217.13,6.837,- 0.0861]}
    A = sch_coeff_dict[gas][0]
    B = sch_coeff_dict[gas][1]
    C = sch_coeff_dict[gas][2]
    D = sch_coeff_dict[gas][3]
    S = A + B*T + C*T**2 + D*T**3
    return S

def k600(Q,V,D,S,eqn=1):
    '''returns the k600 gas exchange velocity for a Schmidt number of 600 (CO2 at 20C) in m/d, given the stream discharge in m3/s, stream velocity
    in m/s, average depth in m, and average slope S.  Uses the equation indexed by the number from table 2 Raymond 2012.'''
    if eqn==1:
        k600 = (V*S)**0.89*D**0.54*5037
        return k600
    if eqn==7:
        k600 = 4725.*(V*S)**0.86*Q**-.14*D**0.6
        return k600
    else:
        print("equation not programmed, no k600 calculated.")
    
def solubility_oil(gas,T,oil_type='light'):
    '''calculates the Henry's coefficient in atm Kg_oil/mol for gas i at temperature T in celcius.  Equations taken from 
    Kharaka and Sprecht 1987 in  Porcelli, D.; Ballentine, C. J. & Wieler, R. (ed.) Noble Gases in Geochemistry and 
    Cosmochemistry Mineralogical Society of America , 2002, 47 pg 547.  Oil is either light or heavy.  Valid for temperatures
    from 25-100C'''
    if oil_type == 'light':    
        sol_dict = {'He':[3.25,-0.0054],
                    'Ne':[3.322,-0.0063],
                    'Ar':[2.121,-0.0003],
                    'Kr':[1.607,0.0019],
                    'Xe':[1.096,0.0035]};
    else:
        sol_dict = {'He':[3.008,-0.0037],
                    'Ne':[2.912,-0.0032],
                    'Ar':[2.03,0.001],
                    'Kr':[1.537,0.0014],
                    'Xe':[0.848,0.0052]};
    A = sol_dict[gas[0:2]][0];
    B = sol_dict[gas[0:2]][1];
    log_K_m = A + B*T;
    K_m = 10.**(log_K_m);
    return K_m


def ng_stm_part_coef(gas,T,S=0):
    """ X_ i = ng_stm_part_coef('gas',T,S,P).  Returns the noble 
    gas steam partion coefficient defined as C_g_i/C_w_i where C_g_i is the 
    concentration of gas_i in the steam phase (mol gas_i/gH2O_vapor) and C_w_i 
    is the concentration of the gas in the water phase (mol gas_i/g_H20_liquid).  
    This partition coeffient is a function of the temperature T (C), and the 
    salinity S (mol/l) at which boiling occurs.""" 
    #calculate ng steam partition coefficient (assume 1atm  dry gas pressure)
    # composition of the dry atmosphere these need to be double checked these are for 4He, 20Ne, 40Ar, 84Kr, 129Xe
    #pp_dict = {'He':5.24E-06,'Ne':1.65E-05,'Ar':9.31E-03,'Kr':6.48E-07,'Xe':2.30E-08};
    #pdb.set_trace();
    z_i = atm_std[gas];
    P_vapor = vapor_pressure(T);
    P_da = 0.000101325;  #1 atm dry gas pressure
    #if P_vapor < 1.01325e-4:
        #print 'Not boiling badness may be occuring....'
    #P_vapor = 10**((.7859+.03477*T_c)/(1+.00412*T_c))*1e-7;  #GPa old way of doing it
    p_i = z_i*1.01325e-4; # GPa - assuming a dry atmosphere at 1 atm
    K_i = solubility(gas,T,S);
    x_i = p_i/K_i; #mol_xe/mol_h2o water
    C_w_i = x_i/18.; #mol_xe/g_h2o;
    C_g_i = z_i*(P_da/P_vapor)*(1./18.); #mol_xe/g_h2o gas
    #calculate the steam phase partition coefficient for xenon
    X_i = C_g_i/C_w_i;          
    return X_i;

#vapor pressure with Antione equation
def vapor_pressure(T):
    '''returns P_vapor the vapor pressure in GPa using the Antione equation, given the temperature T in Celcius.  T must be greater than zero and less than 374 C'''
    if T<=99.0:
        A=8.07131;
        B=1730.63;
        C=233.426;
    else:
        A=8.14019;
        B=1810.94;
        C=244.485;
    P = 10**(A-(B/(C+T))); #mmHG
    P=P/760.*101325;  #Pa
    P_vapor=P/1.0e9;
    return P_vapor
    


#This is the function using freesteam - which is the best, but the library is hard to build and breaks 
#whenever I update.  
#def vapor_pressure(T):
#    """P_vapor = vapor_pressure(T) returns the vapor pressure in GPa for pure water
#    at temperature T in C. Uses IAPWS-IF97 and IAPWS-95 formulation using the freesteam
#    libraries http://freesteam.sourceforge.net/"""
#    T_k = T+273.15; #k
#    P_vapor = psat_T(T_k); #(Pa)
#   P_vapor = P_vapor/1.e9; #in GPa
#   return P_vapor


def lapse_rate(Elev):
    """P_e = lapse_rate(Elev).  returns the atmospheric pressure do elevation Elev (m)
    lapse rate taken from Kip Solomon University of Utah. Pressure in GPa.  
    Assumes a pressure of 1 atm at sea level."""
    P_da = ((1-.0065*Elev/288.15)**5.2561)*0.000101325; #GPa From Kip 
    return P_da

def equil_conc(gas,T,S=0.0,P=0.000101325):
    """ C = equil_conc(gas,T,S,P) returns the atmospheric equilibrium concentration of noble gas (gas) of 
    interest in ccSTP/g.  Gas needs to be in He, Ne, Ar, Kr, Xe.  Temperature (T) in 
    celcius P is atmospheric pressure in GPa, and Salinity is in mol/l.  Pressure can
    be derived from elevation using lapse_rate function"""
    #pp_dict = {'He':5.24E-06,'Ne':1.65E-05,'Ar':9.31E-03,'Kr':6.48E-07,'Xe':2.30E-08};
    #p_i = pp_dict[gas]*P;
    p_i = atm_std[gas]*P;
    K_i = solubility(gas[0:2],T,S);
    x_i = p_i/K_i;  #mol fraction of gas_i in water
    C_i = x_i*(22414./18.);  #ccSTP/g
    return C_i;

def aqueous_conc(gas,T,X,TDG,S=0.0):
    """ C = aqueous_conc(gas,T,X,TDG,S=0.0) returns the aqueous concentration of noble gas (gas) of 
    interest in ccSTP/g, given the gas mixing fraction X.  
    Gas needs to be in He, Ne, Ar, Kr, Xe.  Temperature (T) is water temp  in 
    celcius, TDG is dissolved gas pressure in GPa, and Salinity is in mol/l. """
    K_i = solubility(gas[0:2],T,S);
    P_vapor = vapor_pressure(T)
    P_dry = TDG-P_vapor
    p_i = X*P_dry
    x_i = p_i/K_i;  #mol fraction of gas_i in water
    C_i = x_i*(22414./18.);  #ccSTP/g
    return C_i;

def ce_exc(gas,F,A_e,T,S=0.0,P=0.000101325):
    """ C_ex = ce_exc(gas,F,A_e,T,S,P) returns the concentration of gas (gas) due to closed system
    equilibrium excess air addition. F is the fractionation factor and A_e in the initial volume of 
    entrapped air.  Reference:
    Aeschbach-Hertig, W.; Peeters, F.; Beyerle, U. & Kipfer, R. 
    Paleotemperature reconstruction from noble gases in ground water taking into account 
    equilibration with entrapped air Nature, 2000, 405, 1040-1044"""
    #pp_dict = {'He':5.24E-06,'Ne':1.65E-05,'Ar':9.31E-03,'Kr':6.48E-07,'Xe':2.30E-08};
    #z_i = pp_dict[gas];
    z_i = atm_std[gas];
    C_eq = equil_conc(gas,T,S,P);
    C_ex = ((1-F)*A_e*z_i)/(1+((F*A_e*z_i)/C_eq));
    return C_ex

def rayleigh_frac(gas,T_sep,T_ref,S,Ar_frac):
    """ frac_i = rayleigh_frac(gas,T_sep,T_ref,S,Ar_frac) returns frac_i which the fraction gas/40Ar given as a function
    of the residual argon left in the sample.  Equation taken from 
    Porcelli, D.; Ballentine, C. J. & Wieler, R. (ed.) Noble Gases in Geochemistry and 
    Cosmochemistry Mineralogical Society of America , 2002, 47 pg 554. Where T_ref is the temperature
    of ASW that the sample is referenced too and T_sep is the temperature of vapor phase separation in 
    celcius,S is salinity in mol/l?, and Ar_frac is the residual fraction of argon."""
    K_i = solubility(gas,T_sep,S);
    K_Ar = solubility('Ar',T_sep,S);
    alpha = K_i/K_Ar;
    i_0 = equil_conc(gas,T_ref,S,0.000101325);
    Ar_0 = equil_conc('Ar40',T_ref,S,0.000101325);
    frac_0 = i_0/Ar_0;
    frac_i = frac_0*Ar_frac**(alpha-1);
    return frac_i

def rayleigh_frac_oil_water(gas,T_sep,T_ref,S,Ar_frac,oil_type='light'):
    '''returns the frac_i the is the fraction gas/Ar in the groundwater phase, as a function of the residual
    argon left in the sample.'''
    K_i = solubility(gas,T_sep); #GPa
    K_ar = solubility('Ar',T_sep); #GPa
    K_i_oil_m = solubility_oil(gas,T_sep); #atm Kg/mol
    K_ar_oil_m = solubility_oil('Ar',T_sep); #atm Kg/mol
    K_i_m = K_i*9870./55.6; #taken from porcelli pg 546;
    K_ar_m = K_ar*9780./55.6;
    alpha = (K_i_m*K_ar_oil_m)/(K_ar_m*K_i_oil_m);
    i_0 = equil_conc(gas,T_ref,S,0.000101325);
    Ar_0 = equil_conc('Ar40',T_ref,S,0.000101325);
    frac_0 = i_0/Ar_0;
    frac_i = frac_0*Ar_frac**(alpha-1);
    return frac_i
    
def F_w(C_i,C_Ar,gas,T,S):
    """ F_w_i = F_w(C_i,C_Ar,gas,T,S) returns the fractionation factor for gas_i (C_i/C_Ar)_sample/(C_i/C_Ar)_asw
    where the asw values are the equilibrium solubility concentrations for T the temp in celcius and S 
    the salinity in mol/l?.  Gas concentrations should be in ccSTP/g"""
    C_i_0 = equil_conc(gas,T,S);
    C_Ar_0 = equil_conc('Ar40',T,S);
    R_0 = C_i_0/C_Ar_0;
    R_samp = C_i/C_Ar;
    F_w_i = R_samp/R_0;
    return F_w_i
    
def mix(C_1,C_2,nsteps=1000.):
    """C_mix = mix(C_1,C_2,nsteps=100) returns a vector length nsteps+1 starting with a concentration of 
    C_1 and ending with a concentration of C_2 according to mass balance - C_t = xC_1 + (1-x)C_2 where
    x is the mass fraction of water of concentration C_1"""
    step = (1.-0)/nsteps;
    x = arange(1,step,-step);
    x = append(x,0);
    C_mix = x*C_1 + (1-x)*C_2;
    return C_mix

#this function uses freesteam but is not implemented now because freesteam is so hard to keep going
#def density_h2o(T=20):
#    '''returns the density in kg/m3 of saturated water at temp T assumes pressure is enough to keep all water in
#    liquid phase (for temps above 100)'''
#    s1 = steam_Tx(T+273.15,0.0);
#    rho = s1.rho;
#    return rho

def density_h2o(T=20):
    '''returns the density in kg/m3 of fresh water - taken from kip's spreadsheet - not sure of validity range'''
    A=999.842594;
    B=0.06793952;
    C=-0.00909529;
    D=0.0001001685;
    E=-0.000001120083;
    F=0.000000006536332;
    rho=A+B*T+C*T**2+D*T**3+E*T**4+F*T**5;
    return rho
    
def ccSTP_g2kg_M3(gas,c,T_w = 25.):
    '''converts a noble gas concentration in ccSTP/g_h2o to KgHe/M3_h2o for fresh water at temperature T in C'''
    mol_mas_dict = {'He':4,'Ne':20,'Ar':40,'Ar36':36,'Kr':84,'Xe':132};
    mol_vol = 22414.;
    c_2 = c/mol_vol*mol_mas_dict[gas]*density_h2o(T_w);
    return c_2

def kg_M3_2_ccSTP_g(gas,c,T_w = 25.):
    '''converts a noble gas concentration in ccSTP/g_h2o to KgHe/M3_h2o for fresh water at temperature T in C'''
    mol_mas_dict = {'He':4,'Ne':20,'Ar':40,'Ar36':36,'Kr':84,'Xe':132};
    mol_vol = 22414.;
    c_2 = c*mol_vol/mol_mas_dict[gas]/density_h2o(T_w);
    return c_2

def ccSTP_g2mol_kg(gas,c,T_w=25.):
    '''converts a noble gas concentration in ccSTP/g_h2o to molality (mol/kg_h2o) for fresh water at temperature T (C)'''
    mol_vol = 22414.;
    c_2 = c/mol_vol*1000.;
    return c_2

def activity2mass(decay_const,mol_mass,activity):
    '''Converts an activity in Bq to a mass in grams given the decay constant in seconds and the molar mass of the isotope
    of interest'''
    N = activity/decay_const;
    avos = 6.022e23;
    m = (N/avos)*mol_mass;
    return m
    
def get_UU_samp_id(ifile):
    '''get the sample id for a run, useful making summary'''
    book = xlrd.open_workbook(ifile)
    sheet = book.sheet_by_name('Summary')
    runid = sheet.cell_value(1,1)
    if int(runid[-4:-2])<12:
        samp_id = sheet.cell_value(1,0)
    else:
        samp_id = sheet.cell_value(1,5)
    return samp_id

def import_UU_gas_data(ifile,samp_type='diff_samp'):
    '''data=import_UU_gas_data(ifile) Grab the mol fractions out of a compiled data file from the University of Utah Dissolved Gas Lab.  Good for the data sheet as of Early 2013.  Returns a dictionary containing the sample mole fractions and information.'''
    samp_dict = {}
    mole_frac = {}
    conc = {}
    book = xlrd.open_workbook(ifile)
    sheet = book.sheet_by_name('Summary')
    samp_dict['runid'] = sheet.cell_value(1,1)
    samp_dict['amount'] = sheet.cell_value(1,2)
    if int(samp_dict['runid'][-4:-2])<12:
        samp_dict['samp_id'] = sheet.cell_value(1,0)
    else:
        samp_dict['samp_id'] = sheet.cell_value(1,5)
    print('File Name = ' + ifile)
    print('Sample Id = ' + samp_dict['samp_id'])
    T = float(eval(input('Sampling temperature = ')))
    TDG = float(eval(input('Sample total dissolved gas pressure (mmHg) = ')))
    TDG_GPa = TDG/760*0.000101325

    if int(samp_dict['runid'][-4:-2])<=12:
        if samp_type=='diff_samp':
            names = sheet.row_values(10,start_colx=1,end_colx=24)
            values = sheet.row_values(11,start_colx=1,end_colx=24)
            conc_names=[]
            for i in range(len(names)):
                names[i]=names[i].split()[0]
                conc_names.append(names[i][1::]) 
        elif samp_type=='cu_tube':
            print('cu_tube not yet implemented')
            return
        else:
            print('sample type not available, must be diff_samp or cu_tube')
            return

    if int(samp_dict['runid'][-4:-2])==12:
        if samp_type=='diff_samp':
            names = sheet.row_values(10,start_colx=1,end_colx=27)
            values = sheet.row_values(11,start_colx=1,end_colx=27)
            conc_names=[]
            for i in range(len(names)):
                names[i]=names[i].split()[0]
                conc_names.append(names[i][1::]) 
        elif samp_type=='cu_tube':
            print('sample type not available, must be diff_samp or cu_tube')
            return
        else:
            print('sample type not available, must be diff_samp or cu_tube')
            return
    
    if int(samp_dict['runid'][-4:-2])>=13:
        if samp_type=='diff_samp':
            names = sheet.row_values(11,start_colx=1,end_colx=31)
            values = sheet.row_values(12,start_colx=1,end_colx=31)
            conc_names = sheet.row_values(7,start_colx=1,end_colx=31)
        elif samp_type=='cu_tube':
            names = sheet.row_values(11,start_colx=1,end_colx=31)
            values = sheet.row_values(12,start_colx=1,end_colx=31)
            conc_names = sheet.row_values(7,start_colx=1,end_colx=31)
            conc_values = sheet.row_values(8,start_colx=1,end_colx=31)
        else:
            print('sample type not available, must be diff_samp or cu_tube')
            return
    for i in range(len(names)):
        mole_frac[names[i]]=values[i]
        if samp_type == 'cu_tube':
            conc[conc_names[i]]=conc_values[i]
        if samp_type == 'diff_samp':
            p_vapor = vapor_pressure(T)
            p_i = values[i]*(TDG_GPa-p_vapor)
            gas = names[i][1:3]
            try:
                K_i = solubility(gas,T)
                x_i = p_i/K_i
                conc[conc_names[i]] = x_i*(22414./18.);  #ccSTP/g
            except KeyError:
                conc[conc_names[i]] = 'not available'
    samp_dict['T']=T
    samp_dict['TDG']=TDG
    samp_dict['x_mol']=mole_frac
    samp_dict['conc']=conc
    return samp_dict
            
def make_ng_data_summary():
    '''returns a data_dictionary for a series of samples keyed by the sample id.  Each dictionary element is a sample data dictionary. This is all done with a gui for now, but should be more flexible.'''
    nextone=1
    ddict = {}
    while nextone:
        st = eval(input('Sample Type (1=cu_tube,2=diff_samp) = '))
        if st == 1:
            samp_type='cu_tube'
        else:
            samp_type='diff_samp'
        tkinter.Tk().withdraw()
        ifile = tkinter.filedialog.askopenfilename()
        samp_dict = import_UU_gas_data(ifile,samp_type=samp_type)
        ddict[samp_dict['samp_id']]=samp_dict
        nextone=eval(input('Another File? (1=yes, 0=no) '))

    tkinter.Tk().withdraw()
    outfile = tkinter.filedialog.asksaveasfilename()
    pickle.dump(ddict,open(outfile,'wb'))
    
    return ddict

def RRa(He3,He4):
    '''Calculate the R/Ra ratio given the He3 and He4 mol fraction in the sample. Returns RRa the R/Ra ratio'''
    He3a = get_atm_mol_frac('He3')
    He4a = get_atm_mol_frac('He4')
    RRa = (He3/He4)/(He3a/He4a)
    return RRa
