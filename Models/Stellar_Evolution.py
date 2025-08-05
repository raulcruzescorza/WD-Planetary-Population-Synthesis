import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import re
from astropy import units as u
from astropy.io import ascii
from scipy import integrate
from scipy.integrate import quad
import csv
from astropy.cosmology import Planck18 as cosmo
from astropy.cosmology import z_at_value

# Equation for semimajor axis change under tidal decay and stellar mass loss. See Mustill & Villaver (2012).
# This code was extracted from Sanderson et al. (2022).
def eom(t,x,star,Mpl):
    a = x
    Ms = np.interp(t,star.t,star.Ms)
    Rs = np.interp(t,star.t,star.Rs)
    Menv = np.interp(t,star.t,star.Menv)
    Ls = np.interp(t,star.t,star.Ls)
    mdot = np.interp(t,star.t,star.mdot)
    #mean motion
    n = 2*np.pi*np.sqrt(Ms/(a**3))
    #convective timescale
    tconv = (Menv*Rs*Rs/(star.etaf*Ls))**(1./3.)
        
    freq = (np.pi/(n*star.cf*tconv))**star.gammaf
    if freq > 1:
        f2s = star.fprime
    else:
        f2s = star.fprime*freq
    
    merat = Menv/Ms
    mrat = Mpl/Ms
    
    adot_tide = -merat*(1+mrat)*mrat*(Rs/a)**7*Rs*2*f2s/(9*tconv)
    adot_ml = -a*mdot/(Ms+Mpl)
    
    return adot_tide + adot_ml

class Star:
    
    def __init__(self,t,Ms,Rs,Menv,Ls,mdot,etaf=3,gammaf=2,cf=1,fprime=4.5):
        
        self.t = t
        self.Ms = Ms
        self.Rs = Rs
        self.Menv = Menv
        self.Ls = Ls
        self.mdot = mdot
        self.etaf = etaf
        self.cf = cf
        self.gammaf = gammaf
        self.fprime = fprime

def read_star(file):
    
    Lsol = (1*u.Lsun).decompose().to(u.au**2 * u.Msun / u.yr**3).value #Solar luminosity to code units
    Rsol = (1*u.Rsun).to(u.au).value #Solar radius to code units
    
    data = ascii.read(file,format='csv',names=('Time','Teff','logL','Ms','Rs','Me0','Mee','Md'),delimiter=' ')
    
    Ls = 10**data['logL'] * Lsol

    data.add_column(Ls,name='Ls')
    data['Rs'] = data['Rs']*Rsol
    return data

def inside(t,x,star,mpl):
    return x - np.interp(t,star.t,star.Rs)

inside.terminal = True
inside.direction = -1
