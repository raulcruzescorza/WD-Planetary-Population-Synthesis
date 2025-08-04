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

def ChabrierSalpeter_imf(m):
    A = 0.158/(np.log(10))
    c = np.exp(-(-np.log10(0.08))**2/(2*0.69**2))
    mc = 1.0
    if m < 1.0: # Chabrier
        return A * np.exp(-((np.log10(m)-np.log10(0.08))**2/(2*(0.69)**2)))/m
    else: # Salpeter
        return A * c * m ** -2.35

def Chabrier_IMF(m):
    A = 0.158/(np.log(10))
    mc = 1.0
    
    chabrier = A * np.exp(-((np.log10(m)-np.log10(0.08))**2/(2*(0.69)**2)))/m
    
    return chabrier

def Salpeter_IMF(m):
    A = 0.158/(np.log(10))
    c = np.exp(-(-np.log10(0.08))**2/(2*0.69**2))
    mc = 1.0
    
    salpeter = A * c * m ** -2.35
    
    return salpeter

def IMF_massgen(num_samples):
    
    metallicity = 0.02
    zeta = np.log10(metallicity/0.02)
    alphas = np.array([1.593890E+03, 2.706708E+03,1.466143E+02,4.141960E-02,3.426349E-01,1.949814E+01,4.903830,5.212154E-02,1.312179,8.073972E-01])
    betas = np.array([2.053038E+03,1.483131E+03,-1.048442E+02, 4.564888E-02,0.0, 1.758178,0.0,3.166411E-02,-3.294936E-01,0.0])
    gammas = np.array([1.231226E+03,5.772723E+02,-6.795374E+01,2.958542E-02,0.0,-6.008212,0.0,-2.750074E-03,9.231860E-02,0.0])
    etas = np.array([2.327785E+02,7.411230E+01,-1.391127E+01,5.571483E-03,0.0,-4.470533,0.0,-2.271549E-03,2.610989E-02,0.0])
    coeffs = alphas + betas*zeta + gammas*zeta**2 + etas*zeta**3
    taums = np.array([])
    
    mass_range = np.random.uniform(0.08, 2, num_samples)
    
    
    pdf_values = [ChabrierSalpeter_imf(m) for m in mass_range]
    stellar_masses1 = np.random.choice(mass_range, num_samples, p=pdf_values/np.sum(pdf_values))
    stellar_masses = stellar_masses1[stellar_masses1>0.8]

    
    for i in range(len(stellar_masses)):
        ex = max(0.95, min(0.95 - 0.03*(zeta + 0.30103), 0.99))
        mu = max(0.5, 1.0 - 0.01*max(coeffs[5]/(stellar_masses[i]**coeffs[6]), coeffs[7] + coeffs[8]/(stellar_masses[i]**coeffs[9])))
        tbgb = (coeffs[0] + coeffs[1]*stellar_masses[i]**4 + coeffs[2]*stellar_masses[i]**5.5 + stellar_masses[i]**7)/(coeffs[3]*stellar_masses[i]**2 + coeffs[4]*stellar_masses[i]**7)
        tau = tbgb * max(mu,ex)
        taums = np.append(taums,tau/(10**3))
    
    
    return stellar_masses, taums, stellar_masses1
