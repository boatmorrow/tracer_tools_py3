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

class He_prod_rate:
    '''the production rate for a given rock composition of 3he 4he in given units'''
    def __init__(self,three_He_value=0,four_He_value=0,units='atoms/g_rock/yr',rock=0):
        self.three_He_value = three_He_value;
        self.four_He_value = four_He_value;
        self.units = units;
        self.rock_type = rock;

class rock_type:
    '''this will be a way to pass around different rock typs'''
    def __init__(self,name='none',li_capt_prob=0,total_capt_prob=0,composition=0,porosity=0,density=1.):
        self.name = name;
        self.li_capt_prob = li_capt_prob;
        self.total_capt_prob = total_capt_prob;
        self.composition = composition;
        self.porosity = porosity;
        self.density = density;
        
def define_rock_type(rname='avg_upper_crust',porosity = .3):
    '''Defines a rock_type composition, and production rates.  Production rates are given by the He_prod_rate class
    and are in the default units of atoms/g_rock/yr. A rock type needs a composition of at least 
    U, Th, Li, Na, Mg, Al, Si, C, and Li capture probablility and total capture
    probability.  The capture probability can be calculated from the composition and cross sections in principal.'''
    rock = rock_type(name=rname,porosity=porosity);
    if rname == 'avg_upper_crust':
        #composition mass fractions
        comp_dict = {'Li':2.0e-5,'U':2.8e-6,'Na':2.89e-2,'Mg':1.33e-2,'Al':8.04e-2,'Si':3.09e-1,'C':3.24e-3,'Th':1.07e-5};
        rli_capt_prob = 2.05e-4;  #at some point this could be calculated from the composition
        rtotal_capt_prob = 9.79e-3;
        rdensity = 2.7 # g/cc
        rporosity = 0.1; 
    else:
        print('rock type not available yet')
    rock.composition = comp_dict;
    rock.li_capt_prob = rli_capt_prob;
    rock.total_capt_prob = rtotal_capt_prob;
    rock.porosity = rporosity;
    rock.density = rdensity;
    return rock

def calcHe_prod_rate_whole_rock(rock):
    '''Return the He_prod_rate class for the given rock type rock'''
    HPR = He_prod_rate(rock=rock);
    U = HPR.rock_type.composition['U']*1.e6; #ppm
    Th = HPR.rock_type.composition['Th']*1.e6; #ppm
    Na = HPR.rock_type.composition['Na']*1.e2; #percent
    Mg = HPR.rock_type.composition['Mg']*1.e2;  #percent
    Al = HPR.rock_type.composition['Al']*1.e2;  #percent
    Si = HPR.rock_type.composition['Si']*1.e2;  #percent
    C = HPR.rock_type.composition['C']*1.e2;  #percent
    neut_flux = 0.01 * U * (13.8*Na + 5.4*Mg + 5.0*Al + 1.31*Si + 2.0*C) + 0.01 * Th * \
    (6.0*Na + 2.45*Mg + 2.55*Al + 0.56*Si + 0.83*C) + 0.4788*U; #total neutron flux
    rthree_He_prod = neut_flux*(HPR.rock_type.li_capt_prob/HPR.rock_type.total_capt_prob); #atoms/g_rock/yr
    rfour_He_prod = (3.115e6+1.272e5)*U + 7.710e5*Th;  #atoms/g_rock/yr
    HPR.units = 'atoms/g_rock/yr';
    HPR.three_He_value = rthree_He_prod;
    HPR.four_He_value = rfour_He_prod;
    return HPR

def calcHe_prod_rate_refresh(HPR):
    U = HPR.rock_type.composition['U']*1.e6; #ppm
    Th = HPR.rock_type.composition['Th']*1.e6; #ppm
    Na = HPR.rock_type.composition['Na']*1.e2; #percent
    Mg = HPR.rock_type.composition['Mg']*1.e2;  #percent
    Al = HPR.rock_type.composition['Al']*1.e2;  #percent
    Si = HPR.rock_type.composition['Si']*1.e2;  #percent
    C = HPR.rock_type.composition['C']*1.e2;  #percent
    neut_flux = 0.01 * U * (13.8*Na + 5.4*Mg + 5.0*Al + 1.31*Si + 2.0*C) + 0.01 * Th * \
    (6.0*Na + 2.45*Mg + 2.55*Al + 0.56*Si + 0.83*C) + 0.4788*U; #total neutron flux
    rthree_He_prod = neut_flux*(HPR.rock_type.li_capt_prob/HPR.rock_type.total_capt_prob); #atoms/g_rock/yr
    rfour_He_prod = (3.115e6+1.272e5)*U + 7.710e5*Th;  #atoms/g_rock/yr
    HPR.units = 'atoms/g_rock/yr';
    HPR.three_He_value = rthree_He_prod;
    HPR.four_He_value = rfour_He_prod;

def calc_4He_prod_U_Th(U,Th):
    '''calculates all the production rates for a helium given a U and Th value and assuming average
    crustal composition for everything else.  Returns the Helium_prod_rate class'''
    rock = define_rock_type(rname='avg_upper_crust');
    HPR = calcHe_prod_rate_whole_rock(rock)
    rfour_He_prod = (3.115e6+1.272e5)*U + 7.710e5*Th;  #atoms/g_rock/yr
    HPR.units = 'atoms/g_rock/yr';
    HPR.four_He_value = rfour_He_prod;
    return HPR

def switch_units(HPR,units_desired='atoms/g_rock/yr'):
    '''Convert the units for the helium production rate HPR.  Units need be in string 'He_type/vm/t' and I will update this as we go.  Figure's out what units
    we're in, and what we want, and then switches em.  Assumes complete transfer from rock to water'''
    calcHe_prod_rate_refresh(HPR);
#    pdb.set_trace();
    current_units = HPR.units;
    atoms_to_mols = 1./6.0221415e23;
    atoms_to_ccSTP = 22414./6.0221415e23;
    g_rock_to_cc_rock = 1./HPR.rock_type.density;
    cc_to_m3 = 1./100.**3;
    m3_rock_to_m3_rev = 1/(1.-HPR.rock_type.porosity);
    m3_rev_to_m3_water = HPR.rock_type.porosity;
    mols_to_kg_4He = 4./1000.;
    mols_to_kg_3He = 3./1000.;
    yr_to_s = 3.1556e7;
    atoms_4He_to_kg = atoms_to_mols*mols_to_kg_4He;
    atoms_3He_to_kg = atoms_to_mols*mols_to_kg_3He;
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
                HPR.three_He_value = HPR.three_He_value*atoms_4He_to_kg;
                HPR.four_He_value = HPR.four_He_value*atoms_3He_to_kg;
                a = 'kg_he';
            if He_type_wanted == 'ccSTP':
                #pdb.set_trace();
                HPR.three_He_value = HPR.three_He_value*atoms_to_ccSTP;
                HPR.four_He_value = HPR.four_He_value*atoms_to_ccSTP;
                a = 'ccSTP'
            if He_type_wanted == 'mol':
               #pdb.set_trace();
               HPR.three_He_value = HPR.three_He_value*atoms_to_mols;
               HPR.four_He_value = HPR.four_He_value*atoms_to_mols;
               a = 'mol' 
            
        else:
            print('not implemented yet no converted back to atoms. try again now')
            a = 'atoms';
        
        if vm_current == 'g_rock':
            b = 'g_rock';
            if vm_wanted == 'm3_rev':
                HPR.three_He_value = HPR.three_He_value/g_rock_to_cc_rock/cc_to_m3/m3_rock_to_m3_rev;
                HPR.four_He_value = HPR.four_He_value/g_rock_to_cc_rock/cc_to_m3/m3_rock_to_m3_rev;
                b = 'm3_rev';
            if vm_wanted == 'm3_rock':
                HPR.three_He_value = HPR.three_He_value/g_rock_to_cc_rock/cc_to_m3;
                HPR.four_He_value = HPR.four_He_value/g_rock_to_cc_rock/cc_to_m3;
                b = 'm3_rock';
            if vm_wanted == 'm3_h2o':
                HPR.three_He_value = HPR.three_He_value/g_rock_to_cc_rock/cc_to_m3/m3_rock_to_m3_rev/m3_rev_to_m3_water;
                HPR.four_He_value = HPR.four_He_value/g_rock_to_cc_rock/cc_to_m3/m3_rock_to_m3_rev/m3_rev_to_m3_water;
                b = 'm3_aq';
            if vm_wanted == 'g_h2o':
                HPR.three_He_value = HPR.three_He_value/g_rock_to_cc_rock/m3_rock_to_m3_rev/m3_rev_to_m3_water/cc_water_to_g_water;
                HPR.four_He_value = HPR.four_He_value/g_rock_to_cc_rock/m3_rock_to_m3_rev/m3_rev_to_m3_water/cc_water_to_g_water;
                b = 'g_h2o';
            if vm_wanted == 'kg_h2o':
                HPR.three_He_value = HPR.three_He_value/g_rock_to_cc_rock/m3_rock_to_m3_rev/m3_rev_to_m3_water/cc_water_to_g_water*1000.;
                HPR.four_He_value = HPR.four_He_value/g_rock_to_cc_rock/m3_rock_to_m3_rev/m3_rev_to_m3_water/cc_water_to_g_water*1000.;
                b = 'kg_h2o';
            
        else:
            calcHe_prod_rate_refresh();
            print('not implemented yet no conversion converted back g_rock. try again now')
            b = 'g_rock'
        
        if t_current == 'yr':
            c = 'yr';
            if t_wanted == 's':
                HPR.three_He_value = HPR.three_He_value/yr_to_s;
                HPR.four_He_value = HPR.four_He_value/yr_to_s;
                c = 's';
        else:
            calcHe_prod_rate_refresh();
            print('not implemented yet no conversion converted back yr. try again now')
            c = 'yr';
            
    units_current = a + '/' + b + '/' + c
    print(units_current)
    HPR.units = units_current;

def ccSTP_ccRock_yr2kgHE_m3rev_s(value,porosity):
    #ccSTP to Kg He
    KgHe_ccRock_yr = value/22414*4.0/1000;
    KgHe_m3rev_yr = KgHe_ccRock_yr /( 1/(1-porosity) * 1e-6) ;
    KgHe_m3rev_s = KgHe_m3rev_yr / 3.15567e7;
    return KgHe_m3rev_s

    
