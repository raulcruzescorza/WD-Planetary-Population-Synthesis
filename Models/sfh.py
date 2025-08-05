age_of_the_universe = 13.8
timestep_Gyr = 0.5
Total_number = 5230

def sfh_madau_dickinson(z):
    numerador = 0.015 * (1 + z)**2.7
    denominador = 1 + ((1 + z) / 2.9)**5.6
    return numerador / denominador

#print("Calculating the normaliation factor...")

times_for_integral = np.arange(0, age_of_the_universe, timestep_Gyr)
relative_sfr_per_step = []

for t_age in times_for_integral:
    lookback_time = age_of_the_universe - t_age
    if lookback_time < 0: lookback_time = 0
    try:
        z = z_at_value(cosmo.lookback_time, lookback_time * u.Gyr)
    except Exception:
        z = 1000
    relative_sfr_per_step.append(sfh_madau_dickinson(z))

total_area_under_relative_sfr = np.sum(relative_sfr_per_step)
factor_norm = Total_number / (total_area_under_relative_sfr * timestep_Gyr)
#print(f"Factor de normalización calculado: {factor_norm:.2f}\n")

def sfh_normalizada(age_Gyr):
    """
    Returns the number of stars, per time interval and normalised.
    """
    lookback_time = age_of_the_universe - age_Gyr
    if lookback_time < 0: lookback_time = 0
    try:
        z = z_at_value(cosmo.lookback_time, lookback_time * u.Gyr)
    except Exception:
        z = 1000
        
    # The relative rate multiplied by the normalization factor gives the number of stars per Gyr.
    norm_sfh = sfh_madau_dickinson(z) * factor_norm
    
    # We multiply by the time interval to obtain the number of stars in this interval.
    number_of_stars = norm_sfh * timestep_Gyr
    
    try:
        return int(round(number_of_stars.value))
    except AttributeError:
        # This is just in case the number_of_stars is an integer already.
        return int(round(number_of_stars))

