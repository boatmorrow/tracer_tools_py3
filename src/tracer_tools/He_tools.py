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

class rock_type:
    '''this will be a way to pass around different rock typs'''
    def __init__(self,name='none',li_capt_prob=0,total_capt_prob=0,composition=0,porosity=0,density=1.):
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
            self.elev=0
            self.depth=0 #in g/cm^2
            self.incl=0
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
        atoms_to_mols = 1./6.0221415e23;
        atoms_to_ccSTP = 22414./6.0221415e23;
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
        ELEVATION AND LATTITUDE OF HEIDELBERG RIGHT NOW!'''
        inp_file = self.n_sf_spec
        dfsc = pd.read_csv(inp_file) #read in differential neutron flux at surface for specified location (EXPACS spectra)
        Phi_n = intgrt.trapezoid(dfsc['Neutron'],dfsc['Energy']) #integrate for total flux [cm^-2 s^-1)
        Pno = Phi_n/160*3.1556e7 # 160g/cm^2 attenuation lenght in rock for primary cosmogenic neutrons, per second to per year conversion
        #depth in rock like in Fabryka-Martin
        Lamma_nr = 150 #g/cm2
        K_d = np.exp(-self.APR.depth/Lamma_nr) #exponential decrease with depth
        Pn = K_d*Pno
        self.Pn_ev=Pn # neutrons/g_rock/year
        
    def calcPn_mu(self,K_L=1,K_E=1):
        ''''calculate the neutron flux in n/g_rock/yr for:
                the given rock type at the given depth in rock column (z) in g/cm2. '''
        z = self.APR.depth    
        #calculate stopping rate at depth z
        muk = np.array([[0.8450,1029.6], #relative intensities [-] and attenuation coefficients [g/cm^2]
               [-0.05,161.2],
               [.0205,3000.4]]) # taken from Musy 2023 supplementaries
        
        I_mu0 = 190 #mu/g_rock/yr rate of stopping muons at sea level, high latitudes
        Ikz=0
        #pdb.set_trace()
        for k in range(muk.shape[0]): #depth-dependance based on 3 coefficients
            ikk = muk[k,0]*np.exp(-z/muk[k,1])
            Ikz+=ikk
        
        I_muz = I_mu0*Ikz #depth-dependent stopping rate in #mu/g_rock/yr
        
        #calculate neutron yield per muon
        
        #for fc compound factor
        N_Avog = 6.02214076e23
        MtZt=0
        for e in self.APR.mu_dict:
            Mi = self.composition[e]*N_Avog/self.molM_dict[e] #conversion to atomic concentrations
            Zi = self.Z_dict[e]
            MtZt += Mi*Zi
        
        # Total composition weighted average neutron yield
        Y_n = 0
        for e in self.APR.mu_dict:
            Mi = self.composition[e]*N_Avog/self.molM_dict[e]
            Zi = self.Z_dict[e]
            f_c = Mi*Zi/MtZt
            y_i = self.APR.mu_dict[e][1]
            f_d = self.APR.mu_dict[e][0]
            Y_n += f_c*f_d*y_i
    
        # Pn
        self.Pn_mu = K_L*K_E*I_muz*Y_n #K_L and K_E set to one currently. No correction for altitude and latitude
    
    def  calc_phi_n(self,Lamma_ev=160,Lamma_mu=35):
        '''Calculate the neutron fluxs for primary, muon induced and alpha neutrons in n/cm2/yr, 
            for homogenous neutron production. Using phi_n = Lamma*Pn after Musy 2023.'''
        self.phi_n_ev = Lamma_ev*self.Pn_ev
        self.phi_n_mu = Lamma_mu*self.Pn_mu
        self.phi_n_alpha = Lamma_mu*self.Pn_alpha #should that be lamma_ev? 

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
        lamma_e = 0.581e-10 #partial decay constant of 40K to 40Ar [yr^-1]
        lamma_k = 5.463e-10 # total decay constant of 40K? If yes, should be 5.543e-10 [yr^-1]
        Ma = 39.964 #molar mass of 40Ar [g/mol]
        Na = 6.023e23 #avogadro constant [atoms/mol]
        Xk = 1.176e-4 # 40K/K ratio in the Earth [g_K40/g_K]
        F = Xk*Na/Ma*lamma_e/lamma_k*(np.exp(lamma_k*1)-1) # Production of 40Ar [#atoms] per unit time per gram potassium 
        self.APR.Ar40_value = self.composition['K']*F # Ar40 prodution rate [#atoms/g_rock/yr]
        #pdb.set_trace()
        self.calc39Ar_prod(solar=solar,calc_flux=calc_flux) #function to calculate total Argon production, saves to self.APR.Ar39_value
        self.APR.units ='atoms/g_rock/yr'
    
    def calcP39Ar_n(self,solar='avg',calc_flux=1):
        '''calculate the total 39Ar production from neutrons fora all channels.
        Will calculate neutron production and flux first if desired, for each channel 
        for the given rock type at the depth in g/cm2, elevation (masl) and magnetic
        inclination.'''
        nfile=self.n_sf_spec
        if calc_flux:
            #self.calcPn_ev_JFM(solar=solar)
            self.calcPn_ev() #calling functions to calculate the different neutron production rates
            self.calcPn_mu()
            self.calcPn_alpha()
            self.calc_phi_n() #convert to neutron fluxes (with isotropic assumption & effective attenuation lenghts)
        vflag=1
        #read in spectral data
        inp_file = os.path.join(os.path.dirname(__file__), 'Differential_neutron_flux_normalized.csv')
        dfs = pd.read_csv(inp_file) #neutrons/cm^2/s/MeV
        #read in cosmic spectral flux from excel sheet for heidelberg lattitude and elevation
        inp_file = os.path.join(os.path.dirname(__file__), 'HeidelbergNeutronFlux.csv')
        dfsc = pd.read_csv(inp_file) #neutrons/cm^2/s/MeV
        
        #read in the cross section
        inp_file = os.path.join(os.path.dirname(__file__), 'K39_XS.csv')
        K39xsec = pd.read_csv(inp_file)
        inp_file = os.path.join(os.path.dirname(__file__), 'Ca42_XS.csv')
        Ca42xsec = pd.read_csv(inp_file)
        
        #Potassium production
        
        #truncate to the cross section spectrum
        dfst=dfs[dfs['Energy [MeV]']<K39xsec['Energy [MeV]'].max()] #for K
        dfstc = dfsc[dfsc['Energy']<K39xsec['Energy [MeV]'].max()] #for K
        #calculate normalization (total enegry-integrated fluxes) in units [#n/cm^2/s]
        Phi_n_mu = intgrt.trapezoid(dfs.norm_phi_n_mu_avg,dfs['Energy [MeV]']) #
        Phi_n_alpha = intgrt.trapezoid(dfs.norm_phi_n_alpha_avg,dfs['Energy [MeV]'])
        Phi_n_ev = intgrt.trapezoid(dfsc['Neutron'],dfsc['Energy'])
        #calculate normalized spectrum in units 1/MeV
        phibar_mu = dfst.norm_phi_n_mu_avg/Phi_n_mu #
        #phibar_mu = dfst.norm_phi_n_mu_avg
        phibar_alpha = dfst.norm_phi_n_alpha_avg/Phi_n_alpha
        #phibar_alpha = dfst.norm_phi_n_alpha_avg
        phibar_ev = dfstc['Neutron']/Phi_n_ev
        
        #spectral phi_e for each depth
        try:
            len(self.phi_n_mu) #raises type error if phi_n_mu is an float instead of an array
            # self.phi_n_xy scalar or vector if depth-dependant flux [#n/cm^2/yr],
            # outer product gives array phi_ij with i=depth and j=energy in units [#n/cm^2/yr/MeV]
            phi_e_mu = np.outer(self.phi_n_mu,phibar_mu)
            phi_e_alpha = np.outer(self.phi_n_alpha,phibar_alpha)
            phi_e_ev = np.outer(self.phi_n_ev,phibar_ev)  #
        except TypeError: #if self.phi_n_xy is a scalar
            vflag=0
            phi_e_mu = self.phi_n_mu*phibar_mu
            phi_e_alpha = self.phi_n_mu*phibar_mu
            phi_e_ev = self.phi_n_ev*phibar_ev
            
        #fit a function to xsec
        f = interp1d(K39xsec['Energy [MeV]'],K39xsec['Cross-Section [barns]']*1e-24)
        #interpolate to neutron flux
        sigma_new=f(dfst['Energy [MeV]'])
        sigma_newc=f(dfstc['Energy'])
        #Fold
        #the energy normalized spectral neutron flux
        P39Ar_ev_n = np.zeros_like(self.phi_n_mu) #create empty arrays in the right shape to hold values
        P39Ar_mu_n = np.zeros_like(self.phi_n_mu)
        P39Ar_alpha_n = np.zeros_like(self.phi_n_mu)
        
        #calculate target nuclear concentration.
        Kmfrac = self.composition['K'] #weight fraction g_Ka/g_rock
        K39mfrac = Kmfrac*.932581 #from Reviews Ar dating chapter
        M_K = 39.0983 #molar mass of K g/mol
        Na = 6.02e23 #avogadro constant #atoms/mol
        N_tg = K39mfrac/M_K*Na #atoms K39/g_rock
        if vflag:
            for d in range(phi_e_mu.shape[0]):
                #evaporation
                intgrnd_ev=sigma_newc*phi_e_ev[d,:] #for a particular depth
                I_ev=intgrt.trapezoid(intgrnd_ev,dfstc['Energy']) #calculating the folding integral
                P39Ar_ev_n[d] = N_tg*I_ev # in units #Argon39/g_rock/vr
                #muon
                intgrnd_mu=sigma_new*phi_e_mu[d,:] #for a particular depth
                I_mu=intgrt.trapezoid(intgrnd_mu,dfst['Energy [MeV]'])
                P39Ar_mu_n[d] = N_tg*I_mu
                #alpha
                intgrnd_alpha=sigma_new*phi_e_alpha[d,:] #for a particular depth
                I_alpha=intgrt.trapezoid(intgrnd_alpha,dfst['Energy [MeV]'])
                P39Ar_alpha_n[d] = N_tg*I_alpha
        else:
            intgrnd_ev=sigma_newc*phi_e_ev #for a particular depth
            I_ev=intgrt.trapezoid(intgrnd_ev,dfstc['Energy'])
            P39Ar_ev_n = N_tg*I_ev
            #muon
            intgrnd_mu=sigma_new*phi_e_mu #for a particular depth
            I_mu=intgrt.trapezoid(intgrnd_mu,dfst['Energy [MeV]'])
            P39Ar_mu_n = N_tg*I_mu
            #alpha
            intgrnd_alpha=sigma_new*phi_e_alpha #for a particular depth
            I_alpha=intgrt.trapezoid(intgrnd_alpha,dfst['Energy [MeV]'])
            P39Ar_alpha_n = N_tg*I_alpha
        
        #neglecting ca42 for now.  Should be added some day.
        
        self.APR.P39Ar_alpha_n=P39Ar_alpha_n #saving to the class attributes
        self.APR.P39Ar_ev_n=P39Ar_ev_n
        self.APR.P39Ar_mu_n=P39Ar_mu_n
        self.APR.P39Ar_n = P39Ar_ev_n+P39Ar_mu_n+P39Ar_alpha_n #total prod. rate due to neutrons
        self.APR.units ='atoms/g_rock/yr'
    
    #calculate muon production rate
    def calcP39Ar_mu(self,K_L=1,K_E=1):
        ''''calculate the 39Ar production in atoms/g_rock/yr for:
                the given rock type at the given depth in rock column (z) in g/cm2. '''
        depth = self.APR.depth 
        #calculate stopping rate at depth z
        muk = np.array([[0.8450,1029.6], #relative intensities [-] and attenuation coefficients [g/cm^2]
               [-0.05,161.2],
               [.0205,3000.4]])  # taken from Musy 2023 supplementaries
        
        I_mu0 = 190 #mu/g_rck/yr
        Ikz=0
        #pdb.set_trace()
        for k in range(muk.shape[0]):
            ikk = muk[k,0]*np.exp(-depth/muk[k,1])
            Ikz+=ikk
        
        I_muz = I_mu0*Ikz
        
        #calculate 39Ar yield per muon
        
        #for fc chemical compound factor
        N_Avog = 6.02214076e23
        MtZt=0
        for e in self.APR.mu_dict:
            Mi = self.composition[e]*N_Avog/self.molM_dict[e] #conversion to atomic concentrations
            Zi = self.Z_dict[e]
            MtZt += Mi*Zi
        
        # Total composition weighted average neutron yield
        Y_Ar = 0
        channel_dict = {'K':[.93258,0.015],'Ca':[0.969,0.004]} #isotopic abundace and reation probabbility
        for e in channel_dict.keys():
            Mi = self.composition[e]*N_Avog/self.molM_dict[e] #conversion to atomic concentrations
            Zi = self.Z_dict[e]
            f_c = Mi*Zi/MtZt #chemical compound facotr (for the element)
            f_a = channel_dict[e][0] #abundance of target isotope in element
            f_r = channel_dict[e][1] #reaction probability
            f_d = self.APR.mu_dict[e][0] # nuclear capture probability
            Y_e = f_a*f_c*f_d*f_r #Ar39 yield for channel e [atoms Ar39/stopped muon]
            Y_Ar += Y_e
        
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
        self.APR.Ar39_value=self.APR.P39Ar_n+self.APR.P39Ar_mu
        
    
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
        lamma_e = 0.581e-10 #1/yr partial decay const. for Ar40
        lamma_k = 5.463e-10 #1/yr decay constant of K40
        Ma = 39.964 #molar mass of Ar40
        Na = 6.023e23 #avogadro constant
        Xk = 1.176e-4 #isotope ratio 40K/K
        F = Xk*Na/Ma*lamma_e/lamma_k*(np.exp(lamma_k*age)-1) #as above in Ar40_value except for given age instead of per unit time
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
        t_half = 269.2161339422 #yr-1
        lamma = np.log(2)/t_half #decay constant
        P = self.APR.Ar39_value
        Cseq = P/lamma #depth-dependent if Ar39_value is depth-dependent
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

    
