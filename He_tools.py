# -*- coding: utf-8 -*-
"""
Tools for working with helium

Created on Thu Mar 03 12:16:15 2011

@author: gar279
"""

from noble_gas_tools import *
from numpy import *
#from freesteam import *
import string
import pdb


class rock_type:
    '''this will be a way to pass around different rock typs'''
    def __init__(self,name='none',li_capt_prob=0,total_capt_prob=0,composition=0,porosity=0,density=1.):
        self.name = name;
        self.li_capt_prob = li_capt_prob
        self.total_capt_prob = total_capt_prob
        self.composition = composition
        self.porosity = porosity
        self.density = density
        self.HPR = self.He_prod_rate()
        self.APR = self.Ar_prod_rate()

    class He_prod_rate:
        '''the production rate for a given rock composition of 3he 4he in given units'''
        def __init__(self,three_He_value=0,four_He_value=0,units='atoms/g_rock/yr'):
            self.three_He_value = three_He_value;
            self.four_He_value = four_He_value;
            self.units = units;
    
    class Ar_prod_rate:
        '''The radiogenic production rate for a given rock compoosition.'''
        def __init__(self):
            self.Ar40_value = 0
            self.Ar39_value = 0
            self.Ar40_accum = 0
            self.units = 'atoms/g_rock/yr'

    def switch_units(self,units_desired='atoms/g_rock/yr'):
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
                    a = 'kg_he';
                if He_type_wanted == 'ccSTP':
                    #pdb.set_trace();
                    self.HPR.three_He_value = self.HPR.three_He_value*atoms_to_ccSTP;
                    self.HPR.four_He_value = self.HPR.four_He_value*atoms_to_ccSTP;
                    self.APR.Ar39_value = self.APR.Ar39_value*atoms_to_ccSTP;
                    self.APR.Ar40_value = self.APR.Ar40_value*atoms_to_ccSTP;
                    self.APR.Ar40_accum = self.APR.Ar40_accum*atoms_to_ccSTP;
                    a = 'ccSTP'
                if He_type_wanted == 'mol':
                   #pdb.set_trace();
                   self.HPR.three_He_value = self.HPR.three_He_value*atoms_to_mols;
                   self.HPR.four_He_value = self.HPR.four_He_value*atoms_to_mols;
                   self.APR.Ar39_value = self.APR.Ar39_value*atoms_to_mols;
                   self.APR.Ar40_value = self.APR.Ar40_value*atoms_to_mols;
                   self.APR.Ar40_accum = self.APR.Ar40_accum*atoms_to_mols;
                   a = 'mol' 
                
            else:
                calcHe_prod_rate_whole_rock();
                calcAr_prod_rate_whole_rock()
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
                    b = 'm3_rev';
                if vm_wanted == 'm3_rock':
                    self.HPR.three_He_value = self.HPR.three_He_value/g_rock_to_cc_rock/cc_to_m3;
                    self.HPR.four_He_value = self.HPR.four_He_value/g_rock_to_cc_rock/cc_to_m3;
                    self.APR.Ar39_value = self.APR.Ar39_value/g_rock_to_cc_rock/cc_to_m3;
                    self.APR.Ar40_value = self.APR.Ar40_value/g_rock_to_cc_rock/cc_to_m3;
                    self.APR.Ar40_accum = self.APR.Ar40_accum/g_rock_to_cc_rock/cc_to_m3;
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
                    b = 'g_h2o';
                if vm_wanted == 'kg_h2o':
                    print('warning - assumes all produced gas goes to the water phase')
                    self.HPR.three_He_value = self.HPR.three_He_value/g_rock_to_cc_rock/m3_rock_to_m3_rev/m3_rev_to_m3_water/cc_water_to_g_water*1000.;
                    self.HPR.four_He_value = self.HPR.four_He_value/g_rock_to_cc_rock/m3_rock_to_m3_rev/m3_rev_to_m3_water/cc_water_to_g_water*1000.;
                    self.APR.Ar39_value = self.APR.Ar39_value/g_rock_to_cc_rock/m3_rock_to_m3_rev/m3_rev_to_m3_water/cc_water_to_g_water*1000.;
                    self.APR.Ar40_value = self.APR.Ar40_value/g_rock_to_cc_rock/m3_rock_to_m3_rev/m3_rev_to_m3_water/cc_water_to_g_water*1000.;
                    self.APR.Ar40_accum = self.APR.Ar40_accum/g_rock_to_cc_rock/m3_rock_to_m3_rev/m3_rev_to_m3_water/cc_water_to_g_water*1000.;
                    b = 'kg_h2o';
                
            else:
                calcHe_prod_rate_whole_rock();
                calcAr_prod_rate_whole_rock()
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
                calcHe_prod_rate_whole_rock();
                calcAr_prod_rate_whole_rock()
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
                            'C':3.24e-3,'Th':1.07e-5,'K':2.82e-2} #Ballanine and Burnard Reviews in Geochem, pg 488
            rli_capt_prob = 2.05e-4;  #at some point this could be calculated from the composition
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

    def calcAr_prod_rate_whole_rock(self):
        '''calculate the modern day production rate for a given rock type.'''
        lamma_e = 0.581e-10
        lamma_k = 5.463e-10
        Ma = 39.964
        Na = 6.023e23
        Xk = 1.176e-4
        F = Xk*Na/Ma*lamma_e/lamma_k*(exp(lamma_k*1)-1)
        self.APR.Ar40_value = self.composition['K']*F 
        self.APR.units ='atoms/g_rock/yr'

    def calcAr40_accum(self,age):
        '''calculate the accumulated 40Ar for a given rock for a given age in years.  
            Ar amount and volume can be changed but the age must be in years.'''
        if self.APR.units.split('/')[2] != 'yr':
            print('accumulation only works for age in years')
        lamma_e = 0.581e-10 #1/yr
        lamma_k = 5.463e-10 #1/yr
        lamma = lamma_e+lamma_k
        Ma = 39.964
        Na = 6.023e23
        Xk = 1.176e-4
        F = Xk*Na/Ma*lamma_e/lamma_k*(exp(lamma_k*age)-1)
        self.APR.Ar40_accum = self.composition['K']*F 

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

    
