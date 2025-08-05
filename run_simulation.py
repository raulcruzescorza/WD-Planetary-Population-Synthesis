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


from models.imf import IMF_massgen
from models.planet_generator import generate_planetary_systems
from models.stellar_evolution import Star, eom, read_star, inside
from models.sfh import sfh_normalizada


all_ms_masses = np.array([])
all_ms_lifetimes = np.array([])
all_ms_birth_times = np.array([]) 
all_planetary_systems = {}
has_become_wd = np.array([], dtype=bool)
final_systems = {}
final_all_systems = {}

# Erasing/Creating Test Files
with open("InitialPlanetaryPopulation.txt", "w") as f:
    f.write("--- Initial Exoplanet Population ---\n\n")
with open("FinalPlanetaryPopulation.txt", "w") as f:
    f.write("--- Surviving Planetary Population Around WD ---\n\n")

data_track_1p0 = read_star('Data/agb1p0.dat')
data_track_1p5 = read_star('Data/agb1p5.dat')
data_track_2p0 = read_star('Data/agb2p0.dat')

for t in np.arange(0, age_of_the_universe, timestep_Gyr):
    print(f"--- Simulation Time = {t:.1f} Gyr ---")
    k = len(all_ms_masses)
    num = sfh_normalizada(t)
    
    if num>0:
        new_masses, new_lifetimes, new_masses_g08 = IMF_massgen(num)
        new_planetary_systems = generate_planetary_systems(new_masses)
        all_ms_masses = np.append(all_ms_masses, new_masses)
        all_ms_lifetimes = np.append(all_ms_lifetimes, new_lifetimes)
        all_ms_birth_times = np.append(all_ms_birth_times, [t] * num)
        has_become_wd = np.append(has_become_wd, [False] * num)

        # Añade los nuevos sistemas planetarios al diccionario global usando un índice único
        for i in range(num):
            star_id = k + i
            all_planetary_systems[star_id] = new_planetary_systems[f"Star_{i+1}: {new_masses[i]:.4f} M_Sun"]
    
        # Escribe la información de las estrellas recién nacidas en el archivo inicial
        # Usamos el modo "a" (append) para no borrar el contenido anterior
        with open("Results/InitialPlanetaryPopulation.txt", "a") as f_initial:
            for i in range(num):
                star_id = k + i
                dat = all_planetary_systems[star_id]
                f_initial.write(f"Born at t = {t} Gyr - System {star_id}: M_star = {new_masses[i]:.4f} M_Sun, MS Lifetime: {new_lifetimes[i]:.2f} Gyr\n")
                for j in range(len(dat)):
                    f_initial.write(f"  Planet {j}: Semi-Major axis: {dat[j][0]:.2f} Au, Planet Radius: {dat[j][1]:.2f} R_Earth\n")
                f_initial.write("\n")

    # WD Transformation
    for w in range(len(all_ms_masses)):
        t_birth = all_ms_birth_times[w]
        t_lifetime = all_ms_lifetimes[w]
        
        if t > (t_birth + t_lifetime) and not has_become_wd[w]:
            
            print(f"  > Star {w} (born on t = {t_birth:.1f}, lifetime = {t_lifetime:.2f}) becomes a white dwarf.")
            has_become_wd[w] = True
            
            
            current_star_mass_original = all_ms_masses[w]
            current_star_mass = round(current_star_mass_original, 4)
            
            dat = all_planetary_systems[w]
            
            # "a" (append mode)
            with open("Results/FinalPlanetaryPopulation.txt", "a") as f_final:
                f_final.write(f"--- Stellar System {w} Evolved (M_initial = {current_star_mass} M_Sun) ---\n")
                
                if current_star_mass < 1.25:
                    data = data_track_1p0
                    times = np.linspace(0,1e6,1001)
                elif 1.25 <= current_star_mass < 1.75:
                    data = data_track_1p5
                    times = np.linspace(0,4e6,1001)
                else: 
                    data = data_track_2p0
                    times = np.linspace(0,4e6,1001)
                
                f_final.write(f"  White Dwarf Final Mass: {data['Ms'][-1]:.2f} M_Sun\n")
                
                surviving_planets_list = []
                all_planets_list = []
                
                for j in range(len(dat)):
                    #Planet mass assuming 5510kg/m^3 density.
                    rad = dat[j][1] * 6370e3
                    masa = 4/3 * np.pi * rad**3 * 5510
                    masa1 = masa / 1.99e30 #In solar masses
                    a0 = np.array([dat[j][0]])
                    
                    # Create a "Star" class
                    star = Star(data['Time'], data['Ms'], data['Rs'], data['Me0'], data['Ls'], -data['Md'])
                    
                    # Execute Hannah's code: Tidal+Engulfment Model
                    sol = integrate.solve_ivp(eom, (times[0], times[-1]), a0, method='DOP853', t_eval=times, args=[star, masa1],
                                              rtol=rtol, atol=atol, events=inside)
                    
                    
                    if sol["status"] == 0: # 0 means no engulfment
                        final_a = sol['y'][0][-1]
                        f_final.write(f"  Surviving Planet {j}: a_initial={dat[j][0]:.2f} Au -> a_final={final_a:.2f} Au, Radius={dat[j][1]:.2f} R_Earth\n")
                        #Save data from surviving planets.
                        surviving_planets_list.append((final_a, dat[j][1], dat[j][0])) # (a_final, radio, a_inicial)
                    
                    # All planets saved in all_planets_list
                    all_planets_list.append((dat[j][0], sol['y'][0][-1], dat[j][1], sol["status"]))
                
                final_systems[w] = surviving_planets_list
                final_all_systems[w] = all_planets_list
                
                f_final.write("\n")

print("\n--- Simulation completed ---")
