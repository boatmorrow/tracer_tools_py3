# -*- coding: utf-8 -*-
"""
set of utilities for creating a time indexed pandas dataframe with atmospheric concentrations of environemtnal tracers at a site.  C_t_gas is the anthropogenic gases, e.g. that returned by cfc_tools get_gas_conc, and C_t_isotopes is stable isotopes and tritium (or others) such as that propduced by GNIP.read_gnip from stable_isotope_tools.  The produced data frame can be used to write PFLOTRAN modern recharge constraints using  
"""

import pandas as pd
import tracer_tools.cfc_tools as cfc_tools
import tracer_tools.sf6_tools as sf6_tools
import tracer_tools.noble_gas_tools as noble_gas_tools
import pdb
import datetime
import subprocess


def concat_resampled_time_series(C_t_gas,C_t_isotopes,freq='M',hemisphere='NH'):
    '''returns a time series as a pandas DataRange of all the tracers at the frequency desired. Takes heterogeneous data with different dates and sampling intervals, resamples them.  It then combines them with an outer join, and fills missing data somewhat intelligently. For now inputs are cfc/sf6 pandas data frame and isotopes (stable and tritrium) dataframe indexed by time, but should be a list of the raw data files to process.  return ddff a pandas DataRange for the average concentrations in the atmosphere.  This function is really for making a multi tracer input - used to make a pflotran recharge condition right now.  Might should be moved to pflotran tools...'''
    #the tritium function needs to sussed out better.  I think I like my method for stable isotopes. But could also do the inverse distance weighted method as well.  For now I'm going to move to the isotope method, i.e. just extend in reverse using Vienna.  Probalbly the smartest way to do it is to figure out how to correlate Vienna and others to Bedrichov, then use that function to extend in reverse.  This also would allow for uncertainty analysis.
    #yearly_cfc = pd.date_range(datetime.datetime(C_t_gas.index[0].year,1,1),datetime.datetime(C_t_gas.index[-1].year,1,1),freq=freq)
    #yearly_trit = pd.date_range(datetime.datetime(C_t_isotopes.index[0].year,1,1),datetime.datetime(C_t_isotopes.index[-1].year,1,1),freq=freq)
    
    # the grouped data 
    #cfc_grouped = C_t_gas.groupby(yearly_cfc.asof);
    #trit_grouped = C_t_isotopes.groupby(yearly_trit.asof);
    
    # the mean of the groups
    #df_cfc_sampled = cfc_grouped.mean();
    #df_trit_sampled = trit_grouped.mean();
    df_cfc_sampled = C_t_gas.resample('Y').mean()
    df_trit_sampled = C_t_isotopes.resample('Y').mean()

    #put them together
    ttseries = pd.concat([df_cfc_sampled,df_trit_sampled],axis=1,join='outer')
    #ttseries = pd.concat([df_cfc,df_trit],axis=1,join='outer');
    #now how to deal with the missing data issue!  My thoughts are first padd data series which haven't been updated with the most recent
    #and then fill the data for which we don't have historical data with a value which makes sense
    
    #first pad the data
    ddff = ttseries.fillna(method='pad') #rolls everything forward
    
    #now extend in reverse
    ddff.fillna(method='bfill',inplace=True)
    return ddff

def convert2aqueous(ddff,T,P=1.,S=0.,addHe=True,addAr39=True,addKr81=True,hemisphere='NH'):
    '''converts the atm mixing ratio's in ddff returned from MakeTotalTimeSeries to aqueous concentrations (for gases only)
    as a function of the given temperature pressure and salinity.'''
    ddff_aq = ddff.copy();
    #convert cfc11
    # uses the pandas map method (suuper cool!)
    f = lambda x: cfc_tools.equil_conc_cfc(T,x,cfc=11,P=P,S=S);
    ddff_aq['CFC11'] = ddff['CFC11'+hemisphere].map(f);
    #convert cfc 12
    f = lambda x: cfc_tools.equil_conc_cfc(T,x,cfc=12,P=P,S=S);
    ddff_aq['CFC12'] = ddff['CFC12'+hemisphere].map(f);
    #convert cfc 12
    f = lambda x: cfc_tools.equil_conc_cfc(T,x,cfc=113,P=P,S=S);
    ddff_aq['CFC113'] = ddff['CFC113'+hemisphere].map(f);
    #convert sf6
    f = lambda x: sf6_tools.equil_conc_sf6(T,x,P=P,S=S)
    ddff_aq['SF6'] = ddff['SF6'+hemisphere].map(f)
    
    if addHe:
        ddff_aq['He4'] = noble_gas_tools.equil_conc('He',T);  #ccSTP/gh2o
        ddff_aq['He3'] = noble_gas_tools.equil_conc('He',T)*1.399e-6; #ccSTP/gh2o
    if addAr39:
        ddff_aq['Ar39'] = noble_gas_tools.equil_conc('Ar',T);  #ccSTP/gh2o
        mod_act_ratio=1.78e-6;  #Bq per cc(STP)Ar in atmosphere - from Cook book pg 382
        ar_decay_const=7.06e-6; #day^-1 cook book
        ar_decay_const=ar_decay_const/(24*60*60); #seconds
        mass_ar39 = noble_gas_tools.activity2mass(ar_decay_const,39.0,mod_act_ratio);  #per cc(STP)Ar
        mass_ratio_ar39=mass_ar39*22414./39.95  #mass per mass
        f = lambda x: float(x)*mass_ratio_ar39;        
        ddff_aq['Ar39'] = ddff_aq['Ar39'].map(f);
    if addKr81:
        ddff_aq['Kr81'] = noble_gas_tools.equil_conc('Kr',T) * 1.54e-13
    return ddff_aq

def convert2molality(ddff_aq):
    '''Takes aqueous concentration DataRange and converts it to molality for all the tracers.  This is useful
    for input to models etc...'''
    ddff_aq_molal = ddff_aq.copy();
    f = lambda x: float(x)*1.e-12;
    #convert cfc11
    ddff_aq_molal['CFC11'] = ddff_aq_molal['CFC11'].map(f);
    ddff_aq_molal['CFC12'] = ddff_aq_molal['CFC12'].map(f);
    ddff_aq_molal['CFC113'] = ddff_aq_molal['CFC113'].map(f);
    #convert trit
    ddff_aq_molal['trit TU'] = ddff_aq_molal['trit TU'].map(tu2mol_kg);
    #convert sf6
    f = lambda x: float(x)*1.e-15;
    ddff_aq_molal['SF6'] = ddff_aq_molal['SF6'].map(f);
    #convert He
    f = lambda x: float(x)/22414.*1000.
    ddff_aq_molal['He4']=ddff_aq_molal['He4'].map(f);
    ddff_aq_molal['He3']=ddff_aq_molal['He3'].map(f);
    ddff_aq_molal['Ar39']=ddff_aq_molal['Ar39'].map(f);
    ddff_aq_molal['Kr81']=ddff_aq_molal['Kr81'].map(f)
    #should rename the columns here, so the writing a constraint list works better for plotran
    #this would look something like
    ddff_aq_molal = ddff_aq_molal.rename(columns={'cfc11':'CFC11','cfc12':'CFC12','cfc113':'CFC113','SF6 (ppt)':'SF6','trit TU':'H3'});
    #add helium 3 and 4 to the equation
    return ddff_aq_molal

def tu2mol_kg(x):
    '''converts x tu to molality.'''
    tu2molal = 1./10**18*(2.*6.02e23)/18./6.02e23*1000.;
    molal = float(x)*tu2molal;
    return molal        
