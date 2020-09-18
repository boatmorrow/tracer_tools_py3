# -*- coding: utf-8 -*-
"""
Created on Wed Dec 19 12:06:05 2012

@author: wpgardn
"""
from noble_gas_tools import activity2mass

#calculate the mass ratio 39Ar/Ar
mod_act_ratio=1.78e-6;  #Bq per cc(STP)Ar in atmosphere - from Cook book pg 382
ar_decay_const=7.06e-6; #day^-1 cook book
ar_decay_const=ar_decay_const/(24*60*60);
mass_ar39 = activity2mass(ar_decay_const,39.0,mod_act_ratio);  #per cc(STP)Ar
mass_ratio_ar39=mass_ar39*22414./39.95  #mass per mass