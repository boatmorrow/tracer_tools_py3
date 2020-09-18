# -*- coding: utf-8 -*-
"""
Created on Mon Nov 28 09:07:00 2011

utility function for groundwater hydrology

@author: wpgardn
"""

from iapws import IAPWS97

def perm2hydrocond(k,p=101325,T=20.0):
    '''Calculates the hydraulic conductivity in m/s given the permeability in m2 and the pressure in Pa temperature in celcius.
    Defaults to atmospheric pressure and 20 C'''
    T = 273.15 + T
    S = IAPWS97(T=T,P=p/1.e6) #temperature in Kelvin and pressure in MPa.  This will only work for single phase conditions.  
    mu = S.mu;
    rho = S.rho
    g = 9.80665;
    K = k/(mu/(rho*g))
    return K

def hydrocond2perm(K,p=101325,T=20.0):
    '''Calculates the hydraulic conductivity in m/s given the permeability in m2 and the pressure in Pa temperature in celcius.
    Defaults to atmospheric pressure and 20 C'''
    T = 273.15 + T
    S = IAPWS97(T=T,P=p/1.e6) #temperature in Kelvin and pressure in MPa.  This will only work for single phase conditions.  
    mu = S.mu;
    rho = S.rho
    g = 9.80665;
    k = K*(mu/(rho*g))
    return k

def aperture2perm(b):
    '''k = aperture2perm(b) returns the permeability for a given fracture of aperture b, using the "cubic law" Witherspoon 1980.'''
    k = b**2/12
    return k

def perm2aperture(k):
    '''b = perm2aperture(k) returns the fracture aperture given a permeability using the cubic law Witherspoon 1980.'''
    b = (k*12.)**(1./2)
    return b

def darcy2m2(k):
    '''calculates the permeability in a reasonable unit (m2) given a measurment in Darcy's.  Conversion factor taken
    from Bear.'''
    k = k*9.8697e-13;
    return k;

def eff_frac_cond(K_f,K_m,L,b):
    '''return the effective fracture conductivity (or permeability) for a parallel set of fractures given the fracture conductivity (K_f), matrix conductivity (K_m), fracture spacing (L) and aperture (b).'''
    K_e = (L*K_m+b*K_f)/(L+b)
    return K_e

def eff_frac_por(phi_f,phi_m,L,b):
    '''return the effective fracture porosity for a parallel set of fractures given the fracture porosity (phi_f), matrix porosity (phi_m), fracture spacing (L) and aperture (b).'''
    phi_e = (L*phi_m+b*phi_f)/(L+b)
    return phi_e

def mu_w(p=101325,T=20.0):
    '''calculates the viscosity of water in a single phase system at temp T (C) and pressure P (Pa).'''
    T = 273.15 + T
    S = IAPWS97(T=T,P=p/1.e6)
    mu = S.mu
    return mu

def rho_w(p=101325,T=20.0):
    '''calculates the density of water in a single phase system at temp T (C) and pressure P (Pa).'''
    T = 273.15 + T
    S = IAPWS97(T=T,P=p/1.e6)
    rho = S.rho
    return rho
