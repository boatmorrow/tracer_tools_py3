# -*- coding: utf-8 -*-
"""
Tools for working with helium

Created on Thu Mar 03 12:16:15 2011

@author: gar279
"""

import numpy as np
import tracer_tools.noble_gas_tools as ng
import pdb
import pandas as pd
import scipy.integrate as intgrt
from scipy.interpolate import interp1d
import importlib.resources as pkg_resources
import os
import scipy.constants as scic
import xlwings as xw

#define constants for He_tools

# physical constants 
N_AVOG = scic.N_A

# decay rates [yr^-1]
LAMBDA_AR39 = np.log(2)/268 #NUBASE 2020
LAMBDA_K40_TOTAL = 5.543e-10 #Reviews in Mineralogy and Geochemistry Volume 47, 2002, chapter “K-Ar and Ar-Ar Dating”
LAMBDA_K40_TO_AR40 = 5.808e-11 
LAMBDA_K40_TO_CA40 = 4.962e-10

# from NIST, https://www.nist.gov/pml/atomic-weights-and-isotopic-compositions-relative-atomic-masses (visited on 08/14/2023):
# relative isotopic masses, equivalent to molar mass in g/mol, 
#without uncertainty here
MOL_MASS_AR39 = 38.9643130
MOL_MASS_AR40 = 39.9623831237
MOL_MASS_K39 =  38.9637064864
MOL_MASS_K40 =  39.963998166
MOL_MASS_CA42 = 41.95861783
# standard atomic weights (based on the natural abundance)  
ATOMIC_WEIGHT_AR = 39.948
ATOMIC_WEIGHT_K =  39.0983
ATOMIC_WEIGHT_CA = 40.078 
# isotopic composition ratios (nautral abundance, molar ratios mol/mol)
K39_TO_K_RATIO = 0.932581
K40_TO_K_RATIO = 0.000117
CA42_TO_CA_RATIO = 0.00647


class rock_type:
    '''this will be a way to pass around different rock typs'''
    def __init__(self,name='none',li_capt_prob=0,total_capt_prob=1e-6,composition=0,porosity=0,density=1., excel=1, MC_flag=0 ):
        self.name = name;
        self.li_capt_prob = li_capt_prob
        self.total_capt_prob = total_capt_prob
        self.composition = composition # weight fractions g/g (dimensionless)
        self.porosity = porosity
        self.density = density #g/cm^3
        self.HPR = self.He_prod_rate() 
        self.APR = self.Ar_prod_rate()
        self.Pn_alpha = 0 #neutron production rates for radiogenic (alpha,n) reactions
        self.Pn_ev = 0 # evaporation & spallation (aka primary cosmic ray neutrons)
        self.Pn_mu = 0 # slow muon capture
        self.phi_n_ev = 0 #resulting total neutron fluxes in #n/cm^2/yr
        self.phi_n_mu = 0
        self.phi_n_alpha = 0
        #path to the surface flux spectrum of the current location, default Heidelberg
        self.n_sf_spec = os.path.join(os.path.dirname(__file__), 'HeidelbergNeutronFlux.csv') 
        self.Z_dict = {'Li':3, #proton numbers
                       'U':92,
                       'Na':11,
                       'Mg':12,
                       'Al':13,
                       'Si':14,
                       'C':6,
                       'Th':90,
                       'K':19,
                       'O':8,
                       'F':9,
                       'Ca':20, 
                       'Ti':22,
                       'Fe':26} #Ballentine and Burnard Reviews in Geochem, pg 488, fluorine taken from Sramek 2017
        self.molM_dict = {'O':16.00, #atomic weights in g/mol
                          'H':1.01,
                          'C':12.01,
                          'Na':22.99,
                          'Mg':24.31,
                          'Al':26.98,
                          'Si':28.09,
                          'P':30.97,
                          'K':39.10,
                          'Ca':40.08,
                          'Ti':47.87,
                          'Mn':54.85,
                          'Fe':55.85,
                          'Cl':35.45,
                          'B':10.81,
                          'Sm':150.36,
                          'Gd':157.25,
                          'U':238.03,
                          'Th':232.04,
                          'N':14.01,
                          'Ar':39.95} #taken from Alfimov and Ivy-Ochs 2009, table 2
        self.rel_capt_prob = {'O':1.00, #coulomb-capture probability relative to oxygen  
                          'H':0, #missing value in source 
                          'C':0.43,
                          'Na':1.00,
                          'Mg':0.93,
                          'Al':0.76,
                          'Si':0.84,
                          'P':1.04,
                          'K':1.54,
                          'Ca':1.90,
                          'Ti':2.66,
                          'Mn':2.73,
                          'Fe':3.28,
                          'Cl':1.32,
                          'B':0.25,
                          'Sm':4.40,
                          'Gd':5.80,
                          'U':4.70,
                          'Th':3.00,
                          'N':1.02,
                          'Ar':0} # missing value #taken from Alfimov and Ivy-Ochs 2009, table 2
        # variables to control the excelapplication
        self.excel_file_path = os.path.join(os.path.dirname(__file__), 'EXPACS-eng.xlsx')
        self.excel = excel 
        """Default 1, set to 0 if no excel installation is available on your system. The excel attribute activates the 
        calculation of the surface neutron flux with EXPACS, depending on the other parameters (elev, soilmoisture, w_value, latitude)
        If deactivated the calculation defaults back to Heidelberg neutron spectrum """
        self.MC_flag = MC_flag # Excel application is closed if set to zero
        if self.excel == 1:
            self.app = xw.App()
            self.wb =  self.app.books.open(self.excel_file_path)

    class He_prod_rate:
        '''the production rate for a given rock composition of 3he 4he in given units'''
        def __init__(self,three_He_value=0,four_He_value=0,units='atoms/g_rock/yr'):
            self.three_He_value = three_He_value;
            self.four_He_value = four_He_value;
            self.units = units;
    
    class Ar_prod_rate:
        '''The subsurface production rate for a given rock compoosition.
                Elevation in meters, 
                inclination in degrees positive east,
                depth in rock column - g/cm2, 
        '''
        def __init__(self):
            self.xsec = 'xsec_data/TENDL-2019.csv' #cross section file from TENDL as default
            self.xsec_Ca42 = 'xsec_data/Ca42_XS_TENDL2019.csv'
            # parameters to calculate neutron flux with Expacs
            self.elev=2300 #m
            self.soilmoisture=0.25 #Vol-percent default value
            self.w_value = 50 #Wolf number for solar activity
            self.latitude = 50 #degree latitude positiv = north
 
            #other factors
            self.depth=0 #in g/cm^2
            self.icecover = 0 #g/cm^2
            self.top_shielding= 1 #defaults to no shielding from the topography, multiplied to the total production rate

            self.incl=0 # used for surface flux calculation after JFM 

            #additional attributes used in the uncertainty analysis
            self.K39xsec_noise = np.array([0]) #to be multiplied with the cross-section
            self.Ca42xsec_noise = np.array([0])
            self.stop_rate_noise = 1
            self.Ar_yield_noise = 1
            self.n_yield_noise = 1
            self.P_U_Th_noise = 1
            self.n_mu_shape_noise = np.array([0])
            self.n_alpha_shape_noise = np.array([0])
            

            self.Ar40_value = 0
            self.Ar39_value = 0 #total production value
            self.Ar40_accum = 0
            self.Ar39_accum = 0
            self.P39Ar_alpha_n=0 #Argon 39 production rate values
            self.P39Ar_ev_n=0
            self.P39Ar_mu_n=0
            self.P39Ar_n = 0
            self.P39Ar_mu = 0
            #self.P39Ar_value= 0  redundant, named Ar39_value, see above
            #self.P39Ar_accum = 0
            #set to 1 to include Calcium 42
            self.Ca42 = 1 
            self.Ar39_U_Th_value = 0 #direct calculation with Sramek code, can be used for comparison to P39Ar_alpha_n
            self.Cseq = 0
            self.units = 'atoms/g_rock/yr'
            self.mu_dict = {'C':[.090,0.76], #[fraction of muons adsorbed by nucleus (f_d), average neutron yield per mu (y_n)]
                             'O':[.223,.8],
                             'Na':[.432,1.],
                             'Mg':[.538,.6],
                             'Al':[.582,1.26],
                             'Si':[.671,.086],
                             'K':[.830,1.25],
                             'Ca':[.864,0.75],
                             'Fe':[.906,1.1025]} #taken from Fabryka-Martin 1988, Appendix Table B-2

    def switch_units(self,units_desired='atoms/g_rock/yr'): # self=rock_type object? 
        '''Convert the units for the helium production rate HPR.  Units need be in string 'He_type/vm/t' and I will update this as we go.  Figure's out what units
        we're in, and what we want, and then switches em.  Assumes complete transfer from rock to water'''
        self.calcHe_prod_rate_whole_rock() #switch back to default
        self.calcAr_prod_rate_whole_rock()
        current_units = self.HPR.units
        atoms_to_mols = 1./N_AVOG;
        atoms_to_ccSTP = 22414./N_AVOG;
        g_rock_to_cc_rock = 1./self.density;
        cc_to_m3 = 1./100.**3;
        m3_rock_to_m3_rev = 1/(1.-self.porosity);
        m3_rev_to_m3_water = self.porosity;
        mols_to_kg_4He = 4./1000.;
        mols_to_kg_3He = 3./1000.;
        mols_to_kg_40Ar = 4./1000.;
        mols_to_kg_39Ar = 3./1000.;
        yr_to_s = 3.1556e7;
        atoms_4He_to_kg = atoms_to_mols*mols_to_kg_4He;
        atoms_3He_to_kg = atoms_to_mols*mols_to_kg_3He;
        atoms_40Ar_to_kg = atoms_to_mols*mols_to_kg_40Ar;
        atoms_39Ar_to_kg = atoms_to_mols*mols_to_kg_39Ar;
        cc_water_to_g_water = .998; #at 20C
        if current_units == units_desired:
            print('same units, no conversion')
            return
        else:
            He_type_wanted = units_desired.split('/')[0];
            vm_wanted = units_desired.split('/')[1];
            t_wanted = units_desired.split('/')[2];
            He_type_current = current_units.split('/')[0];
            vm_current = current_units.split('/')[1];
            t_current = current_units.split('/')[2];
            if He_type_current == 'atoms':
                a = 'atoms';
                if He_type_wanted == 'kg_he':
                    #for now means I want kg_he
                    self.HPR.three_He_value = self.HPR.three_He_value*atoms_4He_to_kg;
                    self.HPR.four_He_value = self.HPR.four_He_value*atoms_3He_to_kg;
                    self.APR.Ar40_value = self.APR.Ar40_value*atoms_40Ar_to_kg
                    self.APR.Ar39_value = self.APR.Ar39_value*atoms_39Ar_to_kg
                    self.APR.Ar40_accum = self.APR.Ar40_accum*atoms_40Ar_to_kg
                    self.APR.Ar39_accum = self.APR.Ar39_accum*atoms_39Ar_to_kg
                    a = 'kg_he';
                if He_type_wanted == 'ccSTP':
                    #pdb.set_trace();
                    self.HPR.three_He_value = self.HPR.three_He_value*atoms_to_ccSTP;
                    self.HPR.four_He_value = self.HPR.four_He_value*atoms_to_ccSTP;
                    self.APR.Ar39_value = self.APR.Ar39_value*atoms_to_ccSTP;
                    self.APR.Ar40_value = self.APR.Ar40_value*atoms_to_ccSTP;
                    self.APR.Ar40_accum = self.APR.Ar40_accum*atoms_to_ccSTP;
                    self.APR.Ar39_accum = self.APR.Ar39_accum*atoms_to_ccSTP;
                    a = 'ccSTP'
                if He_type_wanted == 'mol':
                   #pdb.set_trace();
                   self.HPR.three_He_value = self.HPR.three_He_value*atoms_to_mols;
                   self.HPR.four_He_value = self.HPR.four_He_value*atoms_to_mols;
                   self.APR.Ar39_value = self.APR.Ar39_value*atoms_to_mols;
                   self.APR.Ar40_value = self.APR.Ar40_value*atoms_to_mols;
                   self.APR.Ar40_accum = self.APR.Ar40_accum*atoms_to_mols;
                   self.APR.Ar39_accum = self.APR.Ar39_accum*atoms_to_mols;
                   a = 'mol' 
                
            else:
                self.HPR.calcHe_prod_rate_whole_rock();
                self.APR.calcAr_prod_rate_whole_rock()
                print('not implemented yet, converted back to atoms.')
                a = 'atoms';
            
            if vm_current == 'g_rock':
                b = 'g_rock';
                if vm_wanted == 'm3_rev':
                    print('warning - evenly distributes the produced gas between rock and porespace')
                    self.HPR.three_He_value = self.HPR.three_He_value/g_rock_to_cc_rock/cc_to_m3/m3_rock_to_m3_rev;
                    self.HPR.four_He_value = self.HPR.four_He_value/g_rock_to_cc_rock/cc_to_m3/m3_rock_to_m3_rev;
                    self.APR.Ar39_value = self.APR.Ar39_value/g_rock_to_cc_rock/cc_to_m3/m3_rock_to_m3_rev;
                    self.APR.Ar40_value = self.APR.Ar40_value/g_rock_to_cc_rock/cc_to_m3/m3_rock_to_m3_rev;
                    self.APR.Ar40_accum = self.APR.Ar40_accum/g_rock_to_cc_rock/cc_to_m3/m3_rock_to_m3_rev;
                    self.APR.Ar39_accum = self.APR.Ar39_accum/g_rock_to_cc_rock/cc_to_m3/m3_rock_to_m3_rev;
                    b = 'm3_rev';
                if vm_wanted == 'm3_rock':
                    self.HPR.three_He_value = self.HPR.three_He_value/g_rock_to_cc_rock/cc_to_m3;
                    self.HPR.four_He_value = self.HPR.four_He_value/g_rock_to_cc_rock/cc_to_m3;
                    self.APR.Ar39_value = self.APR.Ar39_value/g_rock_to_cc_rock/cc_to_m3;
                    self.APR.Ar40_value = self.APR.Ar40_value/g_rock_to_cc_rock/cc_to_m3;
                    self.APR.Ar40_accum = self.APR.Ar40_accum/g_rock_to_cc_rock/cc_to_m3;
                    self.APR.Ar39_accum = self.APR.Ar39_accum/g_rock_to_cc_rock/cc_to_m3;
                    b = 'm3_rock';
                if vm_wanted == 'm3_h2o':
                    print('warning - assumes all produced gas goes to the water phase')
                    self.HPR.three_He_value = self.HPR.three_He_value/g_rock_to_cc_rock/cc_to_m3/m3_rock_to_m3_rev/m3_rev_to_m3_water;
                    self.HPR.four_He_value = self.HPR.four_He_value/g_rock_to_cc_rock/cc_to_m3/m3_rock_to_m3_rev/m3_rev_to_m3_water;
                    self.APR.Ar39_value = self.APR.Ar39_value/g_rock_to_cc_rock/cc_to_m3/m3_rock_to_m3_rev/m3_rev_to_m3_water;
                    self.APR.Ar40_value = self.APR.Ar40_value/g_rock_to_cc_rock/cc_to_m3/m3_rock_to_m3_rev/m3_rev_to_m3_water;
                    self.APR.Ar40_accum = self.APR.Ar40_accum/g_rock_to_cc_rock/cc_to_m3/m3_rock_to_m3_rev/m3_rev_to_m3_water;
                    b = 'm3_aq';
                if vm_wanted == 'g_h2o':
                    print('warning - assumes all produced gas goes to the water phase')
                    self.HPR.three_He_value = self.HPR.three_He_value/g_rock_to_cc_rock/m3_rock_to_m3_rev/m3_rev_to_m3_water/cc_water_to_g_water;
                    self.HPR.four_He_value = self.HPR.four_He_value/g_rock_to_cc_rock/m3_rock_to_m3_rev/m3_rev_to_m3_water/cc_water_to_g_water;
                    self.APR.Ar39_value = self.APR.Ar39_value/g_rock_to_cc_rock/m3_rock_to_m3_rev/m3_rev_to_m3_water/cc_water_to_g_water;
                    self.APR.Ar40_value = self.APR.Ar40_value/g_rock_to_cc_rock/m3_rock_to_m3_rev/m3_rev_to_m3_water/cc_water_to_g_water;
                    self.APR.Ar40_accum = self.APR.Ar40_accum/g_rock_to_cc_rock/m3_rock_to_m3_rev/m3_rev_to_m3_water/cc_water_to_g_water;
                    self.APR.Ar39_accum = self.APR.Ar39_accum/g_rock_to_cc_rock/m3_rock_to_m3_rev/m3_rev_to_m3_water/cc_water_to_g_water;
                    b = 'g_h2o';
                if vm_wanted == 'kg_h2o':
                    print('warning - assumes all produced gas goes to the water phase')
                    self.HPR.three_He_value = self.HPR.three_He_value/g_rock_to_cc_rock/m3_rock_to_m3_rev/m3_rev_to_m3_water/cc_water_to_g_water*1000.;
                    self.HPR.four_He_value = self.HPR.four_He_value/g_rock_to_cc_rock/m3_rock_to_m3_rev/m3_rev_to_m3_water/cc_water_to_g_water*1000.;
                    self.APR.Ar39_value = self.APR.Ar39_value/g_rock_to_cc_rock/m3_rock_to_m3_rev/m3_rev_to_m3_water/cc_water_to_g_water*1000.;
                    self.APR.Ar40_value = self.APR.Ar40_value/g_rock_to_cc_rock/m3_rock_to_m3_rev/m3_rev_to_m3_water/cc_water_to_g_water*1000.;
                    self.APR.Ar40_accum = self.APR.Ar40_accum/g_rock_to_cc_rock/m3_rock_to_m3_rev/m3_rev_to_m3_water/cc_water_to_g_water*1000.;
                    self.APR.Ar39_accum = self.APR.Ar39_accum/g_rock_to_cc_rock/m3_rock_to_m3_rev/m3_rev_to_m3_water/cc_water_to_g_water*1000.;
                    b = 'kg_h2o';
                
            else:
                self.HPR.calcHe_prod_rate_whole_rock();
                self.APR.calcAr_prod_rate_whole_rock()
                print('not implemented yet no conversion converted back g_rock. try again now')
                b = 'g_rock'
            
            if t_current == 'yr':
                c = 'yr';
                if t_wanted == 's':
                    self.HPR.three_He_value = self.HPR.three_He_value/yr_to_s;
                    self.HPR.four_He_value = self.HPR.four_He_value/yr_to_s;
                    self.APR.Ar39_value = self.APR.Ar39_value/yr_to_s;
                    self.APR.Ar40_value = self.APR.Ar40_value/yr_to_s;
                    c = 's';
            else:
                self.HPR.calcHe_prod_rate_whole_rock();
                self.APR.calcAr_prod_rate_whole_rock()
                print('not implemented yet no conversion converted back yr. try again now')
                c = 'yr';
                
        units_current = a + '/' + b + '/' + c
        print(units_current)
        self.HPR.units = units_current;
        self.APR.units = units_current;
        
    def defined_rock_type(self,rname='avg_upper_crust'):
        '''Defines a rock_type composition in mass fraction, and production rates.  Production rates are given by the He_prod_rate class
        and are in the default units of atoms/g_rock/yr. A rock type needs a composition of at least 
        U, Th, Li, Na, Mg, Al, Si, C, and Li capture probablility and total capture
        probability.  The capture probability can be calculated from the composition and cross sections in principal.'''
        if rname == 'avg_upper_crust':
            #composition mass fractions
            comp_dict = {'Li':2.0e-5,'U':2.8e-6,'Na':2.89e-2,\
                            'Mg':1.33e-2,'Al':8.04e-2,'Si':3.09e-1,\
                            'C':3.24e-3,'Th':1.07e-5,'K':2.82e-2,\
                            'O':4.75e-1,'F':5.57e-4,'Ca':3e-2, \
                            'Ti':3e-3,'Fe':3.15e-2} #Ballentine and Burnard Reviews in Geochem, pg 488, flourine taken from Sramek 2017
            rli_capt_prob = 2.05e-4;  #Hardcoded from Ballentine - at some point this should be calculated from the composition
            rtotal_capt_prob = 9.79e-3;
            rdensity = 2.7 # g/cc
            rporosity = 0.1; 
        if rname == 'VariscanGranite':
            #composition mass fractions - Siebel 1995
            #1. W. Siebel, Constraints on Variscan granite emplacement in north-east Bavaria, Germany: further clues from a petrogenetic study of the Mitterteich granite. Geol. Rundschau 84, 384–398 (1995).
            comp_dict = {'Li':2.5e-4,'U':15e-6,'Na':3.11e-2,\
                            'Mg':0.33e-2,'Al':14.44e-2,'Si':.7241,\
                            'C':3.24e-3,'Th':22e-6,'K':5.06e-2,\
                            'O':4.75e-1,'F':5.57e-4,'Ca':1e-2, \
                            'Ti':3e-3,'Fe':1.9e-2} #Ballentine and Burnard Reviews in Geochem, pg 488, flourine taken from Sramek 2017
            rli_capt_prob = 2.05e-4;  #Hardcoded from Ballentine - at some point this should be calculated from the composition
            rtotal_capt_prob = 9.79e-3;
            rdensity = 2.7 # g/cc
            rporosity = 0.1; 
        else:
            print('rock type not available yet')
        self.composition = comp_dict
        self.li_capt_prob = rli_capt_prob
        self.total_capt_prob = rtotal_capt_prob
        self.porosity = rporosity
        self.density = rdensity
    
    def get_neutron_flux_from_expacs(self):

        # Define the input values which will be passed
        hight = self.APR.elev/1000 #conversion to km
        latitude = self.APR.latitude #degree
        w_value = self.APR.w_value #Wolf number
        soil_water = self.APR.soilmoisture #water volume fraction

        # Define the output range where calculated data will be read 
        output_range = 'D35:E174' # Energy spectrum and differential neutron flux

        # Access the active sheet (you can specify a sheet by name if needed)
        sheet = self.wb.sheets.active
        # Pass input values to the specified range in Excel EXPACS
        sheet['B7'].value = hight #km
        sheet['B8'].value = latitude #degree
        sheet['B10'].value = w_value #solar activity wolf number
        sheet['B14'].value = soil_water #water fraction as a volume fraction

        # Trigger a recalculation
        self.wb.app.calculate()

        expacs_spectrum = sheet[output_range].options(np.array, expand='table').value

        if self.MC_flag == 0:
            self.close_excel()

        # Pass values to Excel and force calculation
        #pass_values_to_expacs_and_calculate(file_path, hight, latitude, w_value, soil_water)

        # Read the calculated data from Excel
        #expacs_spectrum = read_calculated_data_from_expacs(file_path, output_range)

        return expacs_spectrum #returns 2d numpy array with shape  (140,2), float64

    def close_excel(self):
        self.wb.save()
        self.wb.close()
        self.app.quit()

    

    def calcPn_alpha(self):
    
        
        '''Code taken from publication:
            
                O. Šrámek, L. Stevens, W. F. McDonough, S. Mukhopadhyay, and 
                J. Peterson: “Subterranean production of neutrons, 39Ar and 
                21Ne: Rates and uncertainties.” Submitted to Geochim. Cosmochim. Acta
                (arXiv:1509.07436).
            
                Ondřej Šrámek, Charles Univ Prague, ondrej.sramek@gmail.com, Sep 2015
        '''
        wtfrac = np.array([
            self.composition['C'],  #6 carbon
            self.composition['O'],    #8 oxygen
            self.composition['F'],  #9 fluorine
            self.composition['Na'],   #11 sodium
            self.composition['Mg'],   #12 magnesium
            self.composition['Al'],   #13 aluminum
            self.composition['Si'],    #14 silicon
            self.composition['K'],   #19 potassium
            self.composition['Ca'],   #20 calcium
            self.composition['Ti'],  #22 titanium
            self.composition['Fe'],   #26 iron
            self.composition['Th'],  #90 thorium
            self.composition['U']    #92 uranium
            ])
        
        el = np.zeros(15) #array to multiply with elem abund of target nuclies
        el[0] = wtfrac[5] #Al
        el[1] = wtfrac[3] #Na
        el[2] = wtfrac[6] #Si
        el[3] = wtfrac[6] #Si
        el[4] = wtfrac[1] #O
        el[5] = wtfrac[4] #Mg
        el[6] = wtfrac[4] #Mg
        el[7] = wtfrac[2] #F
        el[8] = wtfrac[1] #O
        el[9] = wtfrac[10] #Fe
        el[10] = wtfrac[7] #K
        el[11] = wtfrac[9] #Ti
        el[12] = wtfrac[0] #C
        el[13] = wtfrac[8] #Ca
        el[14] =  1. #SF
        
        pp = np.zeros(3) #array to multiply with elem abund of chain parents
        pp[0] = wtfrac[11] #Th
        pp[1] = wtfrac[12] #U
        pp[2] = wtfrac[12] #U
        
        cn = np.array([
            [2.65E+09, 3.31E+08, 5.03E+09],
            [6.07E+09, 8.02E+08, 1.23E+10],
            [1.95E+08, 2.53E+07, 3.91E+08],
            [1.68E+08, 2.05E+07, 3.17E+08],
            [8.77E+07, 1.33E+07, 2.27E+08],
            [1.72E+09, 2.41E+08, 3.72E+09],
            [1.01E+09, 1.44E+08, 2.22E+09],
            [1.60E+10, 2.34E+09, 3.75E+10],
            [9.50E+06, 1.43E+06, 2.47E+07],
            [1.28E+08, 3.35E+06, 9.55E+07],
            [1.10E+08, 8.92E+06, 1.64E+08],
            [4.35E+08, 2.20E+07, 5.00E+08],
            [3.61E+08, 5.48E+07, 9.96E+08],
            [2.98E+07, 2.22E+06, 4.29E+07],
            [3.01E+03, 2.37E+03, 4.44E+08]
        ])

        nr = np.zeros(cn.shape)
        
        for i in range(cn.shape[0]):
            for j in range(cn.shape[1]):
                nr[i,j] = cn[i,j] * el[i] * pp[j]
        
        Sn = nr.sum() #neutrons/kg_rock/yr
        try:
            len(self.APR.depth)
            self.Pn_alpha = Sn/1000.*np.ones(len(self.APR.depth))
        except TypeError:
            self.Pn_alpha=Sn/1000. #neutrons/g_rock/yr #hand shape issues
        self.Pn_alpha = self.Pn_alpha*self.APR.P_U_Th_noise

    def calcPn_ev_JFM(self,solar='avg',Pno=2000):
        '''Estimates neutron production from high-enery cosmic particles. 
                solar activity one of min, max or average.
                Pno surface production at 70 degrees and 0 masl - 2000n/g_rock/yr JFM
            returns Pn - neutron production in n/g_rock/yr'''
        #geomagnetic latitude lambda_m, latitude factor K_L (solar min), K_L [solar max), atmospheric attenuation length Lambda_n_a
        kl_array = np.array([[0,.164,.233,212],  
                             [10,.170,.241,212],
                             [20, .206,.288,206],
                             [30, .307,.407,195],
                             [40,.511,.617,181],
                             [50,.781,.883,167],
                             [60,.949,1.,164],
                             [70,1.,1.,164],
                             [90,1,1,164]])
        
        #declination effect
        G = self.APR.incl*np.pi/180 # inclination (input attribute) in radiants
        lamma_m = np.arctan(0.5*np.tan(G))*180/np.pi #geomagnetic latitude
        if solar=='min':
            K_L = np.interp(lamma_m,kl_array[:,0],kl_array[:,1])
        if solar=='max':
            K_L = K_L = np.interp(lamma_m,kl_array[:,0],kl_array[:,2])
        else:
            K_L = np.interp(lamma_m,kl_array[:,0],(kl_array[:,1]+kl_array[:,2])/2)
        
        Lamma_na = np.interp(lamma_m,kl_array[:,0],kl_array[:,3])
        
        #elevation effect
        d_atm = ng.lapse_rate(self.APR.elev)*1e9/9.8/10 #g/cm2
        d_sl = ng.lapse_rate(0)*1e9/9.8/10 #g/cm2
        K_E = np.exp(-(d_atm-d_sl)/Lamma_na) 
        
        #depth in rock
        Lamma_nr = 150 #g/cm2
        K_d = np.exp(-self.APR.depth/Lamma_nr)
        Pn = K_E*K_L*K_d*Pno
        self.Pn_ev=Pn #n/g_rock/yr
    
    def calcPn_ev(self):
        '''calculates the production using the excel worksheet spectral flux and equation
        17 from Musy.  inp_file is the spectral energy flux, right now uses
        ELEVATION AND LATTITUDE OF HEIDELBERG if self.APR.excel = 0'''

        if self.excel == 0:
            inp_file = self.n_sf_spec
            dfsc = pd.read_csv(inp_file) #read in differential neutron flux at surface for specified location (EXPACS spectra)
            Phi_n = intgrt.trapezoid(dfsc['Neutron'],dfsc['Energy']) #integrate for total flux [cm^-2 s^-1)
        else:
            dfsc = self.get_neutron_flux_from_expacs()
            Phi_n = intgrt.trapezoid(dfsc[:,1],dfsc[:,0]) #integrate for total flux [cm^-2 s^-1)

        if self.APR.icecover != 0:
            shielding = snow_shielding(self.APR.soilmoisture, self.APR.icecover)
            Phi_n = Phi_n*shielding
        Pno = Phi_n/160*3.1556e7 # 160g/cm^2 attenuation lenght in rock for primary cosmogenic neutrons, per second to per year conversion
        #depth in rock like in Fabryka-Martin
        Lamma_nr = 160 #g/cm2 #150 in JFM, 160 in Gosse-Phillips
        K_d = np.exp(-self.APR.depth/Lamma_nr) #exponential decrease with depth
        Pn = K_d*Pno
        self.Pn_ev=Pn # neutrons/g_rock/year
        
    def calcPn_mu(self,K_L=1,K_E=1):
        ''''calculate the neutron flux in n/g_rock/yr for:
                the given rock type at the given depth in rock column (z) in g/cm2. '''
        z = self.APR.depth + self.APR.icecover    
        #calculate stopping rate at depth z
        muk = np.array([[0.8450,1029.6], #relative intensities [-] and attenuation coefficients [g/cm^2]
               [-0.05,161.2],
               [0.2050,3000.4]]) # taken from Musy 2023 supplementaries
        
        I_mu0 = 190 #mu/g_rock/yr rate of stopping muons at sea level, high latitudes
        Ikz=0
        #pdb.set_trace()
        for k in range(muk.shape[0]): #depth-dependance based on 3 coefficients
            ikk = muk[k,0]*np.exp(-z/muk[k,1])
            Ikz+=ikk
        
        I_muz = I_mu0*Ikz*self.APR.stop_rate_noise #depth-dependent stopping rate in #mu/g_rock/yr
        
        #calculate neutron yield per muon
        
        #for fc compound factor
        MtWt=0
        for e in self.APR.mu_dict:
            Mi = self.composition[e]*N_AVOG/self.molM_dict[e] #conversion to atomic concentrations
            Wi = self.rel_capt_prob[e]
            MtWt += Mi*Wi
        
        # Total composition weighted average neutron yield
        Y_n = 0
        for e in self.APR.mu_dict:
            Mi = self.composition[e]*N_AVOG/self.molM_dict[e]
            Wi = self.rel_capt_prob[e]
            f_c = Mi*Wi/MtWt
            y_i = self.APR.mu_dict[e][1]
            f_d = self.APR.mu_dict[e][0]
            Y_n += f_c*f_d*y_i

        Y_n = Y_n * self.APR.n_yield_noise
        # Pn
        self.Pn_mu = K_L*K_E*I_muz*Y_n #K_L and K_E set to one currently. No correction for altitude and latitude
    
    def  calc_phi_n(self,Lamma_ev=160,Lamma_mu=35):
        '''Calculate the neutron fluxs for primary, muon induced and alpha neutrons in n/cm2/yr, 
            for homogenous neutron production. Using phi_n = Lamma*Pn after Musy 2023.'''
        self.phi_n_ev = Lamma_ev*self.Pn_ev
        self.phi_n_mu = Lamma_mu*self.Pn_mu
        self.phi_n_alpha = Lamma_mu*self.Pn_alpha #lamma_alpha ~ lamma_mu, as neutrons have similar energy range 

    def calcHe_prod_rate_whole_rock(self):
        '''Calculate the He_prod_rate for the rock type '''
        U = self.composition['U']*1.e6; #ppm
        Th = self.composition['Th']*1.e6; #ppm
        Na = self.composition['Na']*1.e2; #percent
        Mg = self.composition['Mg']*1.e2;  #percent
        Al = self.composition['Al']*1.e2;  #percent
        Si = self.composition['Si']*1.e2;  #percent
        C = self.composition['C']*1.e2;  #percent
        neut_flux = 0.01 * U * (13.8*Na + 5.4*Mg + 5.0*Al + 1.31*Si + 2.0*C) + 0.01 * Th * \
        (6.0*Na + 2.45*Mg + 2.55*Al + 0.56*Si + 0.83*C) + 0.4788*U; #total neutron flux
        rthree_He_prod = neut_flux*(self.li_capt_prob/self.total_capt_prob); #atoms/g_rock/yr
        rfour_He_prod = (3.115e6+1.272e5)*U + 7.710e5*Th;  #atoms/g_rock/yr
        self.HPR.units = 'atoms/g_rock/yr';
        self.HPR.three_He_value = rthree_He_prod;
        self.HPR.four_He_value = rfour_He_prod;
        self.HPR.neut_flux = neut_flux

    def calcAr_prod_rate_whole_rock(self,solar='avg',calc_flux=1):
        '''calculate the modern day production rate for a given rock type.  
        The 39Ar production for a given rock type at the given 
        at the elevation (masl), inclination (magnetic), depth (g/cm2), and 
        solar activity (one of (max,min,avg).'''
        lamma_e = LAMBDA_K40_TO_AR40 #partial decay constant of 40K to 40Ar [yr^-1]
        lamma_k = LAMBDA_K40_TOTAL # total decay constant of 40K? If yes, should be 5.543e-10 [yr^-1]
        Ma = MOL_MASS_AR40 #molar mass of 40Ar [g/mol]
        Xk = K40_TO_K_RATIO # 40K/K ratio in the Earth [g_K40/g_K]
        F = Xk*N_AVOG/Ma*lamma_e/lamma_k*(np.exp(lamma_k*1)-1) # Production of 40Ar [#atoms] per unit time per gram potassium 
        self.APR.Ar40_value = self.composition['K']*F # Ar40 prodution rate [#atoms/g_rock/yr]
        #pdb.set_trace()
        self.calc39Ar_prod(solar=solar,calc_flux=calc_flux) #function to calculate total Argon production, saves to self.APR.Ar39_value
        self.APR.units ='atoms/g_rock/yr'
    
    def calcP39Ar_n(self,solar='avg',calc_flux=1):
        '''calculate the total 39Ar production from neutrons fora all channels.
        Will calculate neutron production and flux first if desired, for each channel 
        for the given rock type at the depth in g/cm2, elevation (masl) and magnetic
        inclination.'''
        #nfile=self.n_sf_spec
        if calc_flux:
            #self.calcPn_ev_JFM(solar=solar)
            self.calcPn_ev() #calling functions to calculate the different neutron production rates
            self.calcPn_mu()
            self.calcPn_alpha()
            self.calc_phi_n() #convert to neutron fluxes (with isotropic assumption & effective attenuation lenghts)
        vflag=1 #sets default for depth-dependent calculation
        #read in spectral data: normalized differential flux from Felsenbergkeller (Grieger et al. 2021)
        inp_file = os.path.join(os.path.dirname(__file__), 'Differential_neutron_flux_normalized.csv')
        diff_flux_normal = pd.read_csv(inp_file) #neutrons/cm^2/s/MeV
        # add noise for uncertainty analysis
        if len(diff_flux_normal.norm_phi_n_mu_avg) == len(self.APR.n_mu_shape_noise):
            diff_flux_normal.norm_phi_n_mu_avg += self.APR.n_mu_shape_noise
            #exclude possible negativ values and renew normalisation
            diff_flux_normal.norm_phi_n_mu_avg = diff_flux_normal.norm_phi_n_mu_avg.abs()
            diff_flux_normal.norm_phi_n_mu_avg = diff_flux_normal.norm_phi_n_mu_avg/intgrt.trapezoid(diff_flux_normal.norm_phi_n_mu_avg, diff_flux_normal['Energy [MeV]'])
        if len(diff_flux_normal.norm_phi_n_alpha_avg) == len(self.APR.n_alpha_shape_noise):
            diff_flux_normal.norm_phi_n_alpha_avg += self.APR.n_alpha_shape_noise
            #exclude possible negativ values and renew normalisation
            diff_flux_normal.norm_phi_n_alpha_avg = diff_flux_normal.norm_phi_n_alpha_avg.abs()
            diff_flux_normal.norm_phi_n_alpha_avg = diff_flux_normal.norm_phi_n_alpha_avg/intgrt.trapezoid(diff_flux_normal.norm_phi_n_alpha_avg, diff_flux_normal['Energy [MeV]'])

        #read the input cosmic spectral flux from excel sheet for heidelberg lattitude and elevation (diff_flux_inp)
        if self.excel == 0:
            inp_file = self.n_sf_spec
            diff_flux_inp = pd.read_csv(inp_file) #neutrons/cm^2/s/MeV
        else:
            spectrum_narray = self.get_neutron_flux_from_expacs()
            diff_flux_inp = pd.DataFrame({'Energy':spectrum_narray[:,0], 'Neutron': spectrum_narray[:,1]})

        #read in the cross section
        inp_file = os.path.join(os.path.dirname(__file__), self.APR.xsec)
        K39xsec = pd.read_csv(inp_file, delimiter=';', header=0, names=['energy','xsection'])
        if self.APR.xsec != 'xsec_data/K39_XS_Musy.csv':
            K39xsec.energy = K39xsec.energy*10**-6
            print('Assuming energy in xsec file given in eV. Multiplying e-6 to convert to MeV')
        if len(K39xsec.energy)==len(self.APR.K39xsec_noise):
            K39xsec['xsection'] += self.APR.K39xsec_noise

        inp_file = os.path.join(os.path.dirname(__file__), self.APR.xsec_Ca42)
        Ca42xsec = pd.read_csv(inp_file, delimiter=';',header=0, names=['energy','xsection'])
        if self.APR.xsec_Ca42 != 'xsec_data/Ca42_XS_Musy.csv':
            Ca42xsec.energy = Ca42xsec.energy*10**-6
            print('Assuming energy in Ca42 xsec file given in eV. Multiplying e-6 to convert to MeV')
        if len(Ca42xsec.energy)==len(self.APR.Ca42xsec_noise):
            Ca42xsec['xsection'] +=self.APR.Ca42xsec_noise
        
        #Potassium production
        
        #truncate to the cross section spectrum
        #truncating upper bound
        diff_flux_norm_trunc=diff_flux_normal[diff_flux_normal['Energy [MeV]']<K39xsec['energy'].max()] #for K
        diff_flux_inp_trunc = diff_flux_inp[diff_flux_inp['Energy']<K39xsec['energy'].max()] #for K
        #truncate lower bound
        energy_lowbound = K39xsec['energy'].iloc[0]
        start_index = diff_flux_norm_trunc[diff_flux_norm_trunc['Energy [MeV]']>=energy_lowbound].index[0]
        diff_flux_norm_trunc = diff_flux_norm_trunc.loc[start_index:].reset_index(drop=True)
        start_index = diff_flux_inp_trunc[diff_flux_inp_trunc['Energy']>=energy_lowbound].index[0]
        diff_flux_inp_trunc = diff_flux_inp_trunc.loc[start_index:].reset_index(drop=True)


        #calculate normalization (total enegry-integrated fluxes) in units [#n/cm^2/s]
        Phi_n_mu = intgrt.trapezoid(diff_flux_normal.norm_phi_n_mu_avg,diff_flux_normal['Energy [MeV]']) #
        Phi_n_alpha = intgrt.trapezoid(diff_flux_normal.norm_phi_n_alpha_avg,diff_flux_normal['Energy [MeV]'])
        Phi_n_ev = intgrt.trapezoid(diff_flux_inp['Neutron'],diff_flux_inp['Energy'])
        #calculate normalized spectrum in units 1/MeV (spectral shape)
        phibar_mu = diff_flux_norm_trunc.norm_phi_n_mu_avg/Phi_n_mu #
        #phibar_mu = dfst.norm_phi_n_mu_avg
        phibar_alpha = diff_flux_norm_trunc.norm_phi_n_alpha_avg/Phi_n_alpha
        #phibar_alpha = dfst.norm_phi_n_alpha_avg
        phibar_ev = diff_flux_inp_trunc['Neutron']/Phi_n_ev
        
        #spectral phi_e for each depth
        try:
            len(self.phi_n_mu) #raises type error if phi_n_mu is an float instead of an array
            # self.phi_n_xy scalar or vector if depth-dependant flux [#n/cm^2/yr],
            # outer product gives array phi_ij with i=depth and j=energy in units [#n/cm^2/yr/MeV]
            phi_e_mu = np.outer(self.phi_n_mu,phibar_mu)
            phi_e_alpha = np.outer(self.phi_n_alpha,phibar_alpha)
            phi_e_ev = np.outer(self.phi_n_ev,phibar_ev)  #
        except TypeError: #if self.phi_n_xy is a scalar
            vflag=0 #sets vflag to calculations for a single depth
            phi_e_mu = self.phi_n_mu*phibar_mu
            phi_e_alpha = self.phi_n_mu*phibar_mu
            phi_e_ev = self.phi_n_ev*phibar_ev
            
        #fit a function to xsec
        f = interp1d(K39xsec['energy'],K39xsec['xsection']*1e-24)
        #interpolate to neutron flux
        sigma_new=f(diff_flux_norm_trunc['Energy [MeV]'])
        sigma_newc=f(diff_flux_inp_trunc['Energy'])
        #Fold the energy normalized spectral neutron flux
        P39Ar_ev_n = np.zeros_like(self.phi_n_mu) #create empty arrays in the right shape to hold values
        P39Ar_mu_n = np.zeros_like(self.phi_n_mu)
        P39Ar_alpha_n = np.zeros_like(self.phi_n_mu)
        
        #calculate target nuclear concentration.
        Kmfrac = self.composition['K'] #weight fraction g_Ka/g_rock
        K39mfrac = Kmfrac*.932581 #from Reviews Ar dating chapter (mol fraction)
        M_K = ATOMIC_WEIGHT_K #molar mass of K g/mol
        N_tg = K39mfrac/M_K*N_AVOG #atoms K39/g_rock
        if vflag: #calculating the whole depth range
            for d in range(phi_e_mu.shape[0]): 
                #evaporation
                intgrnd_ev=sigma_newc*phi_e_ev[d,:] #for a particular depth
                I_ev=intgrt.trapezoid(intgrnd_ev,diff_flux_inp_trunc['Energy']) #calculating the folding integral
                P39Ar_ev_n[d] = N_tg*I_ev # in units #Argon39/g_rock/vr
                #muon
                intgrnd_mu=sigma_new*phi_e_mu[d,:] #for a particular depth
                I_mu=intgrt.trapezoid(intgrnd_mu,diff_flux_norm_trunc['Energy [MeV]'])
                P39Ar_mu_n[d] = N_tg*I_mu
                #alpha
                intgrnd_alpha=sigma_new*phi_e_alpha[d,:] #for a particular depth
                I_alpha=intgrt.trapezoid(intgrnd_alpha,diff_flux_norm_trunc['Energy [MeV]'])
                P39Ar_alpha_n[d] = N_tg*I_alpha
        else:
            intgrnd_ev=sigma_newc*phi_e_ev #for a particular depth
            I_ev=intgrt.trapezoid(intgrnd_ev,diff_flux_inp_trunc['Energy'])
            P39Ar_ev_n = N_tg*I_ev
            #muon
            intgrnd_mu=sigma_new*phi_e_mu #for a particular depth
            I_mu=intgrt.trapezoid(intgrnd_mu,diff_flux_norm_trunc['Energy [MeV]'])
            P39Ar_mu_n = N_tg*I_mu
            #alpha
            intgrnd_alpha=sigma_new*phi_e_alpha #for a particular depth
            I_alpha=intgrt.trapezoid(intgrnd_alpha,diff_flux_norm_trunc['Energy [MeV]'])
            P39Ar_alpha_n = N_tg*I_alpha
        
        #Adds the production rates for Ca42(n,alpha)Ar39
        
        if self.APR.Ca42 != 0:
            #print('Including spallation reactions with Ca42.')
            #truncate to the cross section spectrum
            diff_flux_norm_trunc_ca42=diff_flux_normal[diff_flux_normal['Energy [MeV]']<Ca42xsec['energy'].max()] #for Ca
            diff_flux_inp_trunc_ca42 = diff_flux_inp[diff_flux_inp['Energy']<Ca42xsec['energy'].max()] #for Ca
            #truncate lower bound
            energy_lowbound = Ca42xsec['energy'].iloc[0]
            start_index = diff_flux_norm_trunc_ca42[diff_flux_norm_trunc_ca42['Energy [MeV]']>=energy_lowbound].index[0]
            diff_flux_norm_trunc_ca42 = diff_flux_norm_trunc_ca42.loc[start_index:].reset_index(drop=True)
            start_index = diff_flux_inp_trunc_ca42[diff_flux_inp_trunc_ca42['Energy']>=energy_lowbound].index[0]
            diff_flux_inp_trunc_ca42 = diff_flux_inp_trunc_ca42.loc[start_index:].reset_index(drop=True)


            #calculate normalized spectrum in units 1/MeV
            phibar_mu_ca42 = diff_flux_norm_trunc_ca42.norm_phi_n_mu_avg/Phi_n_mu #
            phibar_alpha_ca42 = diff_flux_norm_trunc_ca42.norm_phi_n_alpha_avg/Phi_n_alpha
            phibar_ev_ca42 = diff_flux_inp_trunc_ca42['Neutron']/Phi_n_ev

            #spectral phi_e for each depth for the spectral length of Ca42
            try:
                len(self.phi_n_mu) #raises type error if phi_n_mu is an float instead of an array
                # self.phi_n_xy scalar or vector if depth-dependant flux [#n/cm^2/yr],
                # outer product gives array phi_ij with i=depth and j=energy in units [#n/cm^2/yr/MeV]
                phi_e_mu_ca42 = np.outer(self.phi_n_mu,phibar_mu_ca42)
                phi_e_alpha_ca42 = np.outer(self.phi_n_alpha,phibar_alpha_ca42)
                phi_e_ev_ca42 = np.outer(self.phi_n_ev,phibar_ev_ca42)  #
            except TypeError: #if self.phi_n_xy is a scalar
                vflag=0 #sets vflag to calculations for a single depth
                phi_e_mu_ca42 = self.phi_n_mu*phibar_mu_ca42
                phi_e_alpha_ca42 = self.phi_n_mu*phibar_mu_ca42
                phi_e_ev_ca42 = self.phi_n_ev*phibar_ev_ca42

            #fit a function to xsec
            f_ca42 = interp1d(Ca42xsec['energy'],Ca42xsec['xsection']*1e-24)
            #interpolate to neutron flux
            sigma_new_ca42=f_ca42(diff_flux_norm_trunc_ca42['Energy [MeV]'])
            sigma_newc_ca42=f_ca42(diff_flux_inp_trunc_ca42['Energy'])
            #Fold the energy normalized spectral neutron flux

            #calculate target nuclear concentration.
            Cafrac = self.composition['Ca'] #weight fraction g_Ca/g_rock
            Ca42mfrac = Cafrac*.00647 #currently from Wikipedia, check citable source
            M_Ca = ATOMIC_WEIGHT_CA #molar mass of Ca g/mol 
            N_tg_ca42 = Ca42mfrac/M_Ca*N_AVOG #atoms Ca42/g_rock
            if vflag: #calculating the whole depth range
                for d in range(phi_e_mu_ca42.shape[0]): 
                    #evaporation
                    intgrnd_ev_ca42=sigma_newc_ca42*phi_e_ev_ca42[d,:] #for a particular depth
                    I_ev_ca42=intgrt.trapezoid(intgrnd_ev_ca42,diff_flux_inp_trunc_ca42['Energy']) #calculating the folding integral
                    P39Ar_ev_n[d] += N_tg_ca42*I_ev_ca42 # in units #Argon39/g_rock/vr
                    #muon
                    intgrnd_mu_ca42=sigma_new_ca42*phi_e_mu_ca42[d,:] #for a particular depth
                    I_mu_ca42=intgrt.trapezoid(intgrnd_mu_ca42,diff_flux_norm_trunc_ca42['Energy [MeV]'])
                    P39Ar_mu_n[d] += N_tg_ca42*I_mu_ca42
                    #alpha
                    intgrnd_alpha_ca42=sigma_new_ca42*phi_e_alpha_ca42[d,:] #for a particular depth
                    I_alpha_ca42=intgrt.trapezoid(intgrnd_alpha_ca42,diff_flux_norm_trunc_ca42['Energy [MeV]'])
                    P39Ar_alpha_n[d] += N_tg_ca42*I_alpha_ca42
            else:
                intgrnd_ev_ca42=sigma_newc_ca42*phi_e_ev_ca42 #for a particular depth
                I_ev_ca42=intgrt.trapezoid(intgrnd_ev_ca42,diff_flux_inp_trunc_ca42['Energy'])
                P39Ar_ev_n += N_tg_ca42*I_ev_ca42
                #muon
                intgrnd_mu_ca42=sigma_new_ca42*phi_e_mu_ca42 #for a particular depth
                I_mu_ca42=intgrt.trapezoid(intgrnd_mu_ca42,diff_flux_norm_trunc_ca42['Energy [MeV]'])
                P39Ar_mu_n += N_tg_ca42*I_mu_ca42
                #alpha
                intgrnd_alpha_ca42=sigma_new_ca42*phi_e_alpha_ca42 #for a particular depth
                I_alpha_ca42=intgrt.trapezoid(intgrnd_alpha_ca42,diff_flux_norm_trunc_ca42['Energy [MeV]'])
                P39Ar_alpha_n += N_tg_ca42*I_alpha_ca42 

        self.APR.P39Ar_alpha_n=P39Ar_alpha_n #saving to the class attributes
        self.APR.P39Ar_ev_n=P39Ar_ev_n
        self.APR.P39Ar_mu_n=P39Ar_mu_n
        self.APR.P39Ar_n = P39Ar_ev_n+P39Ar_mu_n+P39Ar_alpha_n #total prod. rate due to neutrons
        self.APR.units ='atoms/g_rock/yr'
    
    #calculate muon production rate
    def calcP39Ar_mu(self,K_L=1,K_E=1):
        ''''calculate the 39Ar production in atoms/g_rock/yr for:
                the given rock type at the given depth in rock column (z) in g/cm2. '''
        depth = self.APR.depth + self.APR.icecover 
        #calculate stopping rate at depth z
        muk = np.array([[0.8450,1029.6], #relative intensities [-] and attenuation coefficients [g/cm^2]
               [-0.05,161.2],
               [0.2050,3000.4]])  # taken from Musy 2023 supplementaries
        
        I_mu0 = 190 #mu/g_rck/yr
        Ikz=0
        #pdb.set_trace()
        for k in range(muk.shape[0]):
            ikk = muk[k,0]*np.exp(-depth/muk[k,1])
            Ikz+=ikk
        
        I_muz = I_mu0*Ikz*self.APR.stop_rate_noise
        
        #calculate 39Ar yield per muon
        
        #for fc chemical compound factor
        MtWt=0
        for e in self.APR.mu_dict:
            Mi = self.composition[e]*N_AVOG/self.molM_dict[e] #conversion to atomic concentrations
            Wi = self.rel_capt_prob[e]
            MtWt += Mi*Wi
        
        # Total composition weighted average neutron yield
        Y_Ar = 0
        channel_dict = {'K':[.93258,0.15],'Ca':[0.969,0.004], 'Ca2':[0.021,0.06]} #isotopic abundace and reation probabbility
        for e in channel_dict.keys():
            if e == 'Ca2':
                el = 'Ca'
            else:
                el = e
            Mi = self.composition[el]*N_AVOG/self.molM_dict[el] #conversion to atomic concentrations
            Wi = self.rel_capt_prob[el]
            f_c = Mi*Wi/MtWt #chemical compound factor (for the element)
            f_a = channel_dict[e][0] #abundance of target isotope in element
            f_r = channel_dict[e][1] #reaction probability
            f_d = self.APR.mu_dict[el][0] # nuclear capture probability
            Y_e = f_a*f_c*f_d*f_r #Ar39 yield for channel e [atoms Ar39/stopped muon]
            Y_Ar += Y_e
        Y_Ar = Y_Ar*self.APR.Ar_yield_noise

        # P39Ar
        self.APR.P39Ar_mu = K_L*K_E*I_muz*Y_Ar #production rate due to direct muon  capture [#Ar39/g_rock/yr]
    
    def calc39Ar_prod(self,solar='avg',calc_flux=1):
        ''' calculate the 39Ar production for a given rock type at the given 
        elevation, inclination, depth, solar activity.  Will only use location and 
        depth to calculate neutron production and flux is calc_flux != 0.
        '''
        #calc the neutron production
        self.calcP39Ar_n(solar=solar,calc_flux=calc_flux)
        #calc the muon production
        self.calcP39Ar_mu()
        self.APR.Ar39_value=(self.APR.P39Ar_n+self.APR.P39Ar_mu)*self.APR.top_shielding
        
    
    def calc39Ar_U_Th_prod(self):
        '''Calculates the 39Ar produces from fissiongenic neutrons in the subsurface. 
           Uses code from:
            O. Šrámek, L. Stevens, W. F. McDonough, S. Mukhopadhyay, and 
            J. Peterson: “Subterranean production of neutrons, 39Ar and 
            21Ne: Rates and uncertainties.” Submitted to Geochim. Cosmochim. Acta
            (arXiv:1509.07436).
            
            Ondřej Šrámek, Charles Univ Prague, ondrej.sramek@gmail.com, Sep 2015 '''
        wtfrac = np.array([
            self.composition['C'],  #6 carbon
            self.composition['O'],    #8 oxygen
            self.composition['F'],  #9 fluorine
            self.composition['Na'],   #11 sodium
            self.composition['Mg'],   #12 magnesium
            self.composition['Al'],   #13 aluminum
            self.composition['Si'],    #14 silicon
            self.composition['K'],   #19 potassium
            self.composition['Ca'],   #20 calcium
            self.composition['Ti'],  #22 titanium
            self.composition['Fe'],   #26 iron
            self.composition['Th'],  #90 thorium
            self.composition['U']    #92 uranium
        ])
        
        el = np.zeros(15) #array to multiply with elem abund of target nuclies
        el[0] = wtfrac[5] #Al
        el[1] = wtfrac[3] #Na
        el[2] = wtfrac[6] #Si
        el[3] = wtfrac[6] #Si
        el[4] = wtfrac[1] #O
        el[5] = wtfrac[4] #Mg
        el[6] = wtfrac[4] #Mg
        el[7] = wtfrac[2] #F
        el[8] = wtfrac[1] #O
        el[9] = wtfrac[10] #Fe
        el[10] = wtfrac[7] #K
        el[11] = wtfrac[9] #Ti
        el[12] = wtfrac[0] #C
        el[13] = wtfrac[8] #Ca
        el[14] =  1. #SF
        
        pp = np.zeros(3) #array to multiply with elem abund of chain parents
        pp[0] = wtfrac[11] #Th
        pp[1] = wtfrac[12] #U
        pp[2] = wtfrac[12] #U
        
        ca = np.array([
            [2.59E+08, 2.36E+07, 3.92E+08],
            [5.04E+08, 4.46E+07, 8.02E+08],
            [3.19E+07, 4.05E+06, 6.46E+07],
            [1.37E+07, 1.35E+06, 2.51E+07],
            [1.97E+07, 2.89E+06, 4.68E+07],
            [3.92E+08, 5.43E+07, 8.13E+08],
            [2.75E+08, 4.04E+07, 6.33E+08],
            [1.78E+09, 2.06E+08, 3.17E+09],
            [2.34E+06, 3.39E+05, 5.70E+06],
            [6.81E+06, 2.53E+04, 1.62E+06],
            [1.06E+07, 4.72E+05, 1.22E+07],
            [5.88E+07, 2.29E+06, 5.89E+07],
            [1.31E+08, 2.10E+07, 4.02E+08],
            [3.96E+06, 2.13E+05, 4.80E+06],
            [2.88E+02, 3.09E+02, 4.73E+07]
        ])

        ar = np.zeros(ca.shape)
        
        for i in range(ca.shape[0]):
            for j in range(ca.shape[1]):
                ar[i,j] = ca[i,j] * el[i] * pp[j] * el[10]
        
        S_Ar = ar.sum() #atoms/kg_rock/yr
        self.APR.Ar39_U_Th_value=S_Ar/1000. #atoms/g_rock/yr
        

    def calcAr40_accum(self,age):
        '''calculate the accumulated 40Ar for a given rock for a given age in years.  
            Ar amount and volume can be changed but the age must be in years.'''
        #really this should always change the units to atoms/g_rock/yr, and then switch back
        unitsi = self.APR.units
        self.switch_units('atoms/g_rock/yr')
        #if self.APR.units.split('/')[2] != 'yr':
        #    print('accumulation only works for age in years')
        #if self.APR.units.split('/')[0] != 'atoms':
        #    print('accumulation only works for amount in atoms, change units to atoms/g_rock/yr')
        lamma_e = LAMBDA_K40_TO_AR40 #1/yr partial decay const. for Ar40
        lamma_k = LAMBDA_K40_TOTAL #1/yr decay constant of K40, magic number from Payton: 5.463e-10, not the same value
        Ma = MOL_MASS_AR40 #molar mass of Ar40
        Xk = K40_TO_K_RATIO #isotope ratio 40K/K
        F = Xk*N_AVOG/Ma*lamma_e/lamma_k*(np.exp(lamma_k*age)-1) #as above in Ar40_value except for given age instead of per unit time
        self.APR.Ar40_accum = self.composition['K']*F 
        self.switch_units(unitsi)

    def calcAr39_accum(self,age,Co=0):
        '''Calculate the accumulated 39Ar for a given rock type for an accumulation age in years.
            Ar accumulated is a determined by the rock's 39Ar production rate (so it should be calculated first),
            and the background concentration (atoms/g_rock/yr) at t = 0.  
            This bg concentration is allowed to be non zero, as the production rate 
            can change...'''
        unitsi = self.APR.units
        self.switch_units('atoms/g_rock/yr')
        #really this should always change the units to atoms/g_rock/yr, and then switch back
        t_half = 269.2161339422 #yr-1 halflife of Argon 39
        #pdb.set_trace()
        lamma = np.log(2)/t_half #conversion to decay constant
        P = self.APR.Ar39_value #total production rate value
        try:
            if len(P) != 1: # production rate is vector-valued, i.e. given depth-dependence
                if len(Co) != 1:
                    if len(age) !=1:                  
                        print('Non-uniform production rate, age, and concentration. Results are not to be trusted')
        except TypeError:
            N = P/lamma+(Co-P/lamma)*np.exp(-lamma*age) # solution for t=age and N(t=0)=Co to ODE dN/dt = -lamma*N + P
    #N = P/lamma+np.outer((Co-P/lamma),np.exp(-lamma*age))
            self.APR.Ar39_accum=N
            self.switch_units(unitsi)
    
    def calc39ArSecEqul(self):
        '''Calculate the secular equilibrium concentration.'''
        unitsi = self.APR.units
        self.switch_units('atoms/g_rock/yr')
        P = self.APR.Ar39_value
        Cseq = P/LAMBDA_AR39 #depth-dependent if Ar39_value is depth-dependent
        self.APR.Cseq = Cseq

def calcHe_prod_rate_refresh(rock_type):
    U = rock_type.composition['U']*1.e6; #ppm
    Th = rock_type.composition['Th']*1.e6; #ppm
    Na = rock_type.composition['Na']*1.e2; #percent
    Mg = rock_type.composition['Mg']*1.e2;  #percent
    Al = rock_type.composition['Al']*1.e2;  #percent
    Si = rock_type.composition['Si']*1.e2;  #percent
    C = rock_type.composition['C']*1.e2;  #percent
    neut_flux = 0.01 * U * (13.8*Na + 5.4*Mg + 5.0*Al + 1.31*Si + 2.0*C) + 0.01 * Th * \
    (6.0*Na + 2.45*Mg + 2.55*Al + 0.56*Si + 0.83*C) + 0.4788*U; #total neutron flux
    rthree_He_prod = neut_flux*(rock_type.li_capt_prob/rock_type.total_capt_prob); #atoms/g_rock/yr
    rfour_He_prod = (3.115e6+1.272e5)*U + 7.710e5*Th;  #atoms/g_rock/yr
    rock_type.HPR.units = 'atoms/g_rock/yr';
    rock_type.HPR.three_He_value = rthree_He_prod;
    rock_type.HPR.four_He_value = rfour_He_prod;


def ccSTP_ccRock_yr2kgHE_m3rev_s(value,porosity):
    #ccSTP to Kg He
    KgHe_ccRock_yr = value/22414*4.0/1000;
    KgHe_m3rev_yr = KgHe_ccRock_yr /( 1/(1-porosity) * 1e-6) ;
    KgHe_m3rev_s = KgHe_m3rev_yr / 3.15567e7;
    return KgHe_m3rev_s

def snow_shielding(soilm,coverh): # based on Jannis Weimar 2022 and Schaller 2002 for muons
    # muon attenuation
    muk = np.array([[0.8450,1030], #relative intensities [-] and attenuation coefficients [g/cm^2]
               [-0.05,160],
               [0.2050,3000]]) # taken from Schaller 2002
    Ikz=0
    for k in range(muk.shape[0]): #depth-dependance based on 3 coefficients
        ikk = muk[k,0]*np.exp(-coverh/muk[k,1])
        Ikz+=ikk

    coverh = coverh*10 #conversion from g/cm^2 to mm SWE to use Jannis Weimar's empirical formula
    # soilmoisture correction
    a0 = 3.48 #+-0.48
    a1 = 0.21 #+-0.02  
    N0 = (a0 + a1*soilm)/(a0+soilm)
    # soilmoisture dependence of attenuation of diffusive flux
    a2= 20.2 #+-1.1 mm
    a3=10.9 #+-0.3 mm
    a4=0.12 #+-0.1 
    lam_n = a2 - (a2-a3)*np.exp(-a4*soilm) #mm
    # fractions
    Nm = 0.019 #+-0.001
    Nn = 0.468 #+-0.007
    Nh = 0.513 #+- 0.001
    #attenuation of high-energy flux
    lamh = 1341 #+-15 mm
    # total shielding factor
    Nf = N0*(Nm*Ikz + Nh*np.exp(-coverh/lamh) + Nn*np.exp(-coverh**0.7/lam_n)) 
    return Nf


# defining the functions to access the expacs file and calculate the neutron flux
def pass_values_to_expacs_and_calculate(file_path, hight, latitude, w_value, soil_water):
    # Connect to the Excel application
    app = xw.App()

    # Open the workbook
    wb = app.books.open(file_path)

    # Access the active sheet (you can specify a sheet by name if needed)
    sheet = wb.sheets.active

    # Pass input values to the specified range in Excel EXPACS
    sheet['B7'].value = hight #km
    sheet['B8'].value = latitude #degree
    sheet['B10'].value = w_value #solar activity wolf number
    sheet['B14'].value = soil_water #water fraction as a volume fraction

    # Trigger a recalculation
    wb.app.calculate()

    # Save the workbook (optional, depending on your needs)
    wb.save()

    # Close the workbook
    wb.close()

    # Quit the Excel application
    app.quit()

def read_calculated_data_from_expacs(file_path, output_range):
    # Connect to the Excel application
    app = xw.App()

    # Open the workbook
    wb = app.books.open(file_path)

    # Access the active sheet (you can specify a sheet by name if needed)
    sheet = wb.sheets.active

    # Read the calculated data from the specified range in Excel
    calculated_data = sheet[output_range].options(np.array, expand='table').value

    # Close the workbook
    wb.close()

    # Quit the Excel application
    app.quit()

    return calculated_data
