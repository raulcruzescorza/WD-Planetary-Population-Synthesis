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

def parse_range(s):
    """
    Given a chain of the form "0.50-1.00", return a tuple (min, max).
    """
    s = s.strip()
    partes = s.split('-')
    if len(partes) == 2:
        try:
            return float(partes[0]), float(partes[1])
        except ValueError:
            return None, None
    return None, None
def parse_efficiency_asymmetric(s):
    """
    Processes an efficiency chain and returns the value and the upper and lower error limits.
    
    Returns:
        tuple: (value, upper_error, bottom_error)
    """
    s = s.strip()
    # LaTeX Format: ${2.5}_{-1.7}^{+2.6}\times {10}^{-3}$
    pattern_latex = r'^\$\{([\d\.]+)\}_\{-([\d\.]+)\}\^\{\+([\d\.]+)\}\\times\s*\{10\}\^\{(-?\d+)\}\$$'
    m = re.match(pattern_latex, s)
    if m:
        val, err_neg, err_pos, exponent = map(float, m.groups())
        factor = 10 ** exponent
        value = val * factor
        # Returns the error limits separately
        positive_error = err_pos * factor
        negative_error = err_neg * factor
        return value, positive_error, negative_error

    # Limit format: <5.3 x 10^-4
    pattern_limit = r'^<\s*([\d\.]+)\s*x\s*10\^(-?\d+)$'
    m = re.match(pattern_limit, s)
    if m:
        val, exponent = map(float, m.groups())
        value = val * (10 ** exponent)
        # Para límites superiores, no hay un error definido, devolvemos None
        return value, None, None

    return None, None, None

def leer_tabla_y_procesar(file_path):
    """
    Read an ASCII file and extracts:
      - The original columns "Period", "Radius" y "Combined Detection and Vetting Efficiency"
      - From "Period" it generates "Period_min" and "Period_max"
      - From the "Radius" it generates "Radius_min" and "Radius_max"
      - From the frequency column it extracts "Eff_value" and "Eff_error"
    
    It is assumed that the data lines are separated by tabulations or by two or more spaces.
    """
    data = []
    with open(file_path, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            # Ignore empty lines, headings or comments.
            if not line or line.startswith("Table") or line.startswith("Note"):
                continue
            # Ignore the headings with keywords.
            if re.search(r'Period\s+', line) and re.search(r'Radius\s+', line):
                continue
            # Separate the lines using tabulators or 2 or more spaces.
            columnas = re.split(r'\t+|\s{2,}', line)
            if len(columnas) < 3:
                continue
            # Extract the first three columns: Period, Radius, y Combined frequency
            period_str = columnas[0]
            radius_str = columnas[1]
            eff_str = columnas[2]
            data.append([period_str, radius_str, eff_str])
    
    # Create an initial DataFrame
    df = pd.DataFrame(data, columns=["Period", "Radius", "Combined Efficiency"])
    
    # Process the period column: separate in min and max
    period_min, period_max = [], []
    for s in df["Period"]:
        pmin, pmax = parse_range(s)
        period_min.append(pmin)
        period_max.append(pmax)
    df["Period_min"] = period_min
    df["Period_max"] = period_max
    
    # Process the radius column: separate in min and max
    radius_min, radius_max = [], []
    for s in df["Radius"]:
        rmin, rmax = parse_range(s)
        radius_min.append(rmin)
        radius_max.append(rmax)
    df["Radius_min"] = radius_min
    df["Radius_max"] = radius_max


    eff_value, eff_error_pos, eff_error_neg = [], [], []
    for s in df["Combined Efficiency"]:
        val, err_p, err_n = parse_efficiency_asymmetric(s)
        eff_value.append(val)
        eff_error_pos.append(err_p)
        eff_error_neg.append(err_n)
        
    df["Eff_value"] = eff_value
    df["Eff_error_pos"] = eff_error_pos # New column for the upper error.
    df["Eff_error_neg"] = eff_error_neg # New column for the lower error.

    return df


def generate_planetary_systems(stellar_masses):
    df = leer_tabla_y_procesar("/Users/Raul/Desktop/ajab31abt2_ascii.txt")
    systems = {}
    for i in range(len(stellar_masses)):
        m_star = stellar_masses[i]
        system = []
        for index, row in df.iterrows():
            eff_value = row['Eff_value']
            err_pos = row['Eff_error_pos']
            err_neg = row['Eff_error_neg']
            
            if pd.isna(err_pos) or pd.isna(err_neg):
                current_rate = eff_value
            else:
                # Calculamos la sigma como el promedio de los errores
                sigma = (err_pos + err_neg) / 2.0
                # Tomamos una muestra de una distribución normal
                current_rate = np.random.normal(loc=eff_value, scale=sigma)

            # Nos aseguramos de que la tasa nunca sea negativa
            if current_rate < 0:
                current_rate = 0
            
            # Usamos esta nueva tasa (variable) para la distribución de Poisson
            num = np.random.poisson(current_rate)

            for j in range(num):
                period = np.random.uniform(float(row['Period_min']), float(row['Period_max']))
                a = period_to_distance(period, m_star)
                radius = np.random.uniform(float(row['Radius_min']), float(row['Radius_max']))
                system.append((a, radius))

        systems[f"Star_{i+1}: {m_star:.4f} M_Sun"] = system
        
    return systems
