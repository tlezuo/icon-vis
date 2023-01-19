import math
import numpy as np
from matplotlib import cm

## CLASS DEFINITIN ##
# 1) variable dictionary
class Variable_dictionary:
    def __init__(self, name, title, units, min, max, ticks, cmap, nlev, vemin, vemax, veticks):
        self.name = name
        self.title = title
        self.units = units
        self.min = min
        self.max = max
        self.ticks = ticks
        self.cmap=cmap
        self.modlev=nlev
        self.vemin = vemin
        self.vemax = vemax
        self.veticks = veticks
    
    # name in icon output
    name = 'NAME'
    # title in plot
    title = 'TITLE'
    # color bar extent
    min = 0,100
    max = 100
    ticks = 10
    # colormap plt
    cmap = 'binary'

####################################################################
##  SET DICTS ##
# 2D vars
# 2m T
T_2M = Variable_dictionary('T_2M','2 m temperature','[°C]', 5, 31, 5,'RdYlBu',0,0,0,0)
# 2m Q
QV_2M = Variable_dictionary('qv_2m','2 m specific humidity','[g/kg]', 0,0.016, 0.001, 'YlGnBu_r',0,0,0,0)
# 10m wind
U_10M = Variable_dictionary('u_10m','10m zonal wind','[m/s]',-20,20,2,'RdPu',80,-30,30,2)
V_10M = Variable_dictionary('v_10m','10m meridional wind','[m/s]',-20,20,2,'RdPu',80,-30,30,2)
VEL_10M = Variable_dictionary('VEL_10M','10 m wind speed','[m/s]',0,8.5,1,'RdPu',80,0,41,1)
DIR_10M = Variable_dictionary('DIR_10M','10 m wind  direction','[°]',0,361,10,'viridis',80,0,361,10)

# tke
TKEs = Variable_dictionary('TKEs','Turbulent Kinetic Energy','[J]', 0,1,0.1, 'spring',81,0,1000,50)

# htopo
HSURF = Variable_dictionary('HSURF','topography','[m asl]',0,4000,100,'terrain',0,0,0,0)

# 3D vars
# Altitude
ALT = Variable_dictionary('ALT', 'Altitude', '[m AMSL]', 200, 5000, 250, 'RdYlBu',80,0,0,0)

# Temperature
T = Variable_dictionary('T', 'Temperature', '[°K]', 285,305.5,1, 'YlOrRd',80,270, 311, 1) #levels for vertical

# Pressure
P = Variable_dictionary('P', 'Pressure', '[hPa]', 900, 1000, 10, 'RdYlBu',80,300,1000,50)

# Wind
U = Variable_dictionary('U','Zonal wind speed','[m/s]',-20,20,2,'RdPu',80,-30,30,2)
V = Variable_dictionary('V','Meridional wind speed','[m/s]',-20,20,2,'RdPu',80,-30,30,2)
W = Variable_dictionary('W','Vertical wind speed','[m/s]',-2,2,0.25,'coolwarm',81,-3,3,0.25)
VEL = Variable_dictionary('VEL','Wind speed','[m/s]',0,15,1,'RdPu',80,0,41,1)
DIR = Variable_dictionary('DIR','Wind direction','[°]',0,361,10,'viridis',80,0,361,10)

# Humidity
QV = Variable_dictionary('QV','Specific Humidity','[kg/kg]', 0,0.016, 0.0005, 'YlGnBu',80,0,0.016, 0.0005)
RH = Variable_dictionary('RH','Relative Humidity','[%]', 0,101,5, 'Blues',80,0,101,5,)

# TKE
TKE = Variable_dictionary('TKE','Turbulent Kinetic Energy','[J]', 0,1000,50, 'spring',81,0,1000,50)

# calc vars
# Potential Temperature
TH = Variable_dictionary('TH', 'Potential Temperature', '[°K]', 270, 295, 5, 'YlOrRd',80,270, 295, 5)
TD = Variable_dictionary('TD', 'Dewpoint Temperature', '[°K]', 270, 295, 5, 'YlOrRd',80,270, 295, 5)

# diff vars
TH_diff = Variable_dictionary('TH_diff', 'Potential Temperature change (hourly normalized)', '[°K/h]', -7, 7, 0.5, 'YlOrRd',80,270, 295, 5)
T_diff = Variable_dictionary('T_diff', 'Temperature change (hourly normalized)', '[°K/h]', -7, 7, 0.5, 'YlOrRd',80,270, 295, 5)
VEL_diff = Variable_dictionary('VEL_diff', 'Wind speed change (hourly normalized)', '[°K/h]', -7, 7, 0.5, 'YlOrRd',80,270, 295, 5)

# gradients
TH_grad = Variable_dictionary('TH_grad', 'Potential T gradient (Munich-Kolsass)', '[°K/100km]', -6, 6, 0.5, 'seismic',80,270, 295, 5)
T_grad = Variable_dictionary('T_grad', 'T gradient (Munich-Kolsass)', '[°K/100km]', -6, 6, 0.5, 'seismic',80,270, 295, 5)
P_grad = Variable_dictionary('P_grad', 'Pressure gradient (Munich-Kolsass)', '[hPa/100km]', -2, 2, 0.2, 'BrBG',80,270, 295, 5)
VEL_grad = Variable_dictionary('VEL_grad', 'Velocity gradient (Munich-Kolsass)', '[ms/100km]', -6, 6, 0.5, 'PRGn',80,270, 295, 5)

# Radiation
LW_d = Variable_dictionary('lw_down', 'Incoming longwave radiation', '[W/m2]', 0, 200, 0.5, 'YlOrRd',80,270, 295, 5)
LW_u = Variable_dictionary('lw_up', 'Outgoing longwave radiation', '[W/m2]', 0, 200, 0.5, 'YlOrRd',80,270, 295, 5)
SW_d = Variable_dictionary('sw_down', 'Incoming shortwave radiation', '[W/m2]', 0, 200, 0.5, 'YlOrRd',80,270, 295, 5)
SW_u = Variable_dictionary('sw_up', 'Outgoing shortwave radiation', '[W/m2]', 0, 200, 0.5, 'YlOrRd',80,270, 295, 5)
####################################################################
## CALC FUNCTIONS ##
def calculate_wind_vel_from_uv(u, v): # from pp
    """Calculate wind velocity from U, V components.

    Args:
        u (pd series) u wind component in m/s
        v (pd series) v wind component in m/s

    Returns:
        pd series: wind velocity in m/s

    """
    wind_vel = np.sqrt(u**2 + v**2)

    return wind_vel

def calculate_wind_dir_from_uv(u, v, modulo_180=False): # from pp
    """Calculate wind direction from U, V components.

    Args:
        u (pd series):     u wind component in m/s
        v (pd series):     v wind component in m/s
        modulo_180 (bool): if True, retruned angle will be between [-180,180]

    Returns:
        pd series: wind direction in °

    """
    
    # convert to wind direction coordinate, different from trig unit circle coords
    # if the wind directin is 360 then returns zero (by %360)
    # inspired from wind_uv_to_dir function in:
    # https://github.com/blaylockbk/Ute_WRF/blob/master/functions/wind_calcs.py

    wind_dir = (270 - np.rad2deg(np.arctan2(v, u))) % 360

    # if requested convert from [0,360] to [-180,180]
    if modulo_180 == True:
        wind_dir = (wind_dir + 180) % 360 - 180
        # wind_dir  = np.rad2deg(np.arctan2(v, u))

    return wind_dir

def calculate_potT(T, P):
    """Calculate potential temperature from T and P.

    Args:
        T (pd series):     temperature in K
        P (pd series):     pressure in hPa
        modulo_180 (bool): if True, retruned angle will be between [-180,180]

    Returns:
        pd series: wind direction in °

    """
    #defining parameters
    press_r = 1000                  #reference pressure (hPa)
    Rd = 287                        #specific gas constant for dry air (J/kg*K)
    cp = 1004                       #speific heat of dry air at constant pressure (J/kg*K)

    #compute potential temperature
    potT = T*(press_r/P)**(Rd/cp)

    return potT


def calculate_tdew_from_rh(rh, T, temperature_metric="celsius", verbose=False):
    """Calculate dew point temperature from relative humidity and temperature.

    Args:
        rh (pd series):    air relative humidity in %
        T (pd series):     air temperature in °C
        temperature_metric (str, optional): Input temperature unit. Defaults to "celsius".

    Returns:
        pandas series: dew point temperature timeseries (in °C or K)

    """
    if verbose:
        print(
            "Calculating dew point temperature (dewp_temp) from relative humidity and temp."
        )

    if temperature_metric != "celsius":
        if verbose:
            print("Assuming input temperature unit is Kelvin for rh calculation.")
        T = T - 273  # K to °C

    # inspired from humidity.to.dewpoint in:
    # https://github.com/geanders/weathermetrics/blob/master/R/moisture_conversions.R
    Tdew = (rh / 100) ** (1 / 8) * (112 + (0.9 * T)) - 112 + (0.1 * T)  # in °C

    if temperature_metric != "celsius":
        Tdew = Tdew + 273  # °C to K

    return Tdew

def calculate_qv_from_tdew(Press, Tdew, verbose=False):
    """Calculate specific humidity from pressure and dew point temperature.

    Args:
        Press (pd series):   air pressure in hPa
        Tdew (pd series):    dew point temperature in °C

    Returns:
        pandas series: specific humidity series in kg/kg

    """
    if verbose:
        print("Calculating specific humidity (qv) from press and dewp_temp.")

    # after eq. 4.24 in Practical Meteorology from Stull
    # P in hPa, Td in °C and qv in kg/kg
    e = 6.112 * np.exp((17.67 * Tdew) / (Tdew + 243.5))
    qv = (0.622 * e) / (Press - (0.378 * e))

    return qv

def calculate_qv_from_rh(Press, rh, T, verbose=False):
    """Calculate specific humidity from pressure, relative humidity and temperature.

    Args:
        Press (pd series):   air pressure series in hPa
        rh    (pd series):   air relative humidity in %
        T     (pd series):   air temperature in K

    Returns:
        pandas series: specific humidity series in kg/kg

    """
    if verbose:
        print("Calculating specific humidity (qv) from press and relative humidity.")

    Tdew = calculate_tdew_from_rh(rh, T)

    qv = calculate_qv_from_tdew(Press, Tdew)

    return qv
