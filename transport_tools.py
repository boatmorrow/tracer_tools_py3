''' functionality for transport calculations'''
import numpy as N

def calc_characterstic_diffusion_length(t,D_aq,n,tort):
    '''calculate the characteristic diffusion length l for a given time.  Calculates De as D_aq*n*tort.'''
    De = D_aq*n*tort
    l = N.sqrt(4*De*t)
    return l
