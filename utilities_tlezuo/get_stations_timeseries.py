"""
Functions to read in data from weather stations of ACINN, Innsbruck iBox stations

Author: Tobia Lezuo
Based on: TS_read_stations.py in wrf_crossinn_toni package by Antonia Fritz

Created: 27.10.2022
"""

import numpy as np
import pandas as pd
import datetime
import os.path as path
import glob2
from ipdb import set_trace

def yy2yyyy(yy):
    """Add '20' to timestamp.

    E.g. '21081812' will become '2021081812'

    Args:
        yy (string)

    Returns:
        yyyy (string)

    """
    
    return f"20{yy}"

def parse_timestamps(timestamps):
    """Parse different formats of timestamps.

    Args:
    timestamps  (str, datetime-obj or list of one of those)

    Returns:
    t1, t2      (YYYYmmddHH)

    """
    # if list
    if isinstance(timestamps, list):

        # only 1 element in list
        if len(timestamps) == 1:
            ts0 = timestamps[0]

            # if string
            if isinstance(ts0, str):
                if len(ts0) == 10:
                    return ts0, ts0
                elif len(ts0) == 8:
                    ts0_20 = yy2yyyy(ts0)
                    return ts0_20, ts0_20

            # if datetime object
            elif isinstance(ts0, dt.datetime):
                ts0_str = ts0.strftime("%Y%m%d%H")
                return ts0_str, ts0_str

        # 2 elements in list
        elif len(timestamps) == 2:
            ts0 = timestamps[0]
            ts1 = timestamps[1]
            # if string
            if isinstance(ts0, str):
                if len(ts0) == 10:
                    return ts0, ts1
                elif len(ts0) == 8:
                    ts0_20 = yy2yyyy(ts0)
                    ts1_20 = yy2yyyy(ts1)
                    return ts0_20, ts1_20

            # if datetime-obj
            elif isinstance(ts0, dt.datetime):
                ts0_str = ts0.strftime("%Y%m%d%H")
                ts1_str = ts1.strftime("%Y%m%d%H")
                return ts0_str, ts1_str

        else:
            print("! too many timestamps")
            sys.exit(1)

    # if only 1 string is given
    elif isinstance(timestamps, str):
        if len(timestamps) == 10:
            return timestamps, timestamps
        elif len(timestamps) == 8:
            ts_20 = yy2yyyy(timestamps)
            return ts_20, ts_20

    # if only 1 datetime object is given
    elif isinstance(timestamps, dt.datetime):
        ts_str = timestamps.strftime("%Y%m%d%H")
        return ts_str, ts_str

    else:
        print("! timestamps input is nonsense!")
        sys.exit(1)

# %% Function to parse_timestamps again, adding minutes
def parse_timestamps_min(t):
    """
    Put timestamps to format str ('yyyy-mm-dd HH:MM:SS')

    E.g. '2021081812' will become '2021-08-18 12:00:00'

    Args:
        t10 (string)

    Returns:
        t17 (string)

    """

    t17=str(t[0:4]) +'-'+ str(t[4:6]) +'-'+ str(t[6:8]) +' '+ str(t[8:10]) +':00:00'
    
    return t17

def getKey(dct,value): # finds 
    '''
    Function to find the dict key = station specific variable name for a secific var = value

    Parameters
    ----------
    dct : dictionary
        "var_s":"var"
    value : str
        value=var
        
    Returns
    -------
    var_stationname : str
        key=var_s
    '''
    var_stationname = str([key for key in dct if (dct[key] == value)])
    var_stationname = var_stationname[2:-2] # cuts the edges of the string
    return var_stationname

# %% Function to read in ACINN data
def read_acinn(loc, vars, path_ACINN, start_time,end_time,
               correct_direction=False):
    '''
    Function to read in the station data of ACINN

    Parameters
    ----------
    station : str
        Abreviation for the station as indicated in stations.py in sdf
    vars : list
        short name of desired variable
    path_ACINN : str
        Path to the ACINN data
    timestamps : list of strings     
        either 1 or 2 timestamps YYYYmmddHH
        then converted to start_time, end_time

    Returns
    -------
    data : pandas Dataframe
        Daraframe containing the ACINN data for:
        - temp:    2m Temperature in K
        - press:  surface pressure in hPa
        - wind_vel:  wind speed in m/s
        - wind_dir: wind direction in degree
        - h_m:   height of the wind sensor in m
        The measurements closest to the WRF output heights (agl) are selected.

    '''
    
    # parse timestamp input to 2 string of format YYYYmmddHH
    # t1,t2 = parse_timestamps(timestamps)
    # # parse timestamp to match csv row name of format 'yyyy-mm-dd HH:MM:SS'
    # start_time = parse_timestamps_min(t1)
    # end_time = parse_timestamps_min(t2)

    # As the ACINN data sets aren't labelled consistently, a dict is needed
    # to bring them all to the same names. Only important parameters are kept.
    # Description of parameters can be found at:
    # https://acinn-data.uibk.ac.at/pages/station-list.html
    # https://acinn-data.uibk.ac.at/pages/station-list.html
    acinn_v = {'hoch': {
                        # 'ta_avg': 'T_2M',   # RAW: Air temperature, ventilation [°C] 
                        'tair2':'T_2M', #	air temperature (PT100) 1.5m?
                        # 'TAIR2I':'T_2M' #air temperature at sonic 2 7.08m
                        'p': 'P', # RAW: Air pressure [hPa]
                        'meanu2':'VEL_10M', # FLUXL: mean rotated and unfiltered u wind component (streamwise) [m/s] sonic2 7.08m
                                            # is VEL as streamwise MEANV=0
                        # 'w_dir_avg':'DIR_10M', # RAW: wind dir [deg] at 1m
                        'wind_dir2':'DIR_10M', #FLUXL: wind direction [deg] at 7.08m
                        'cnr4_lw_in_wm2':'lw_down', # RADIAT1
                        'cnr4_lw_out_wm2':'lw_up', # RADIAT1
                        'cnr4_sw_in_wm2':'sw_down', # RADIAT1
                        'cnr4_sw_out_wm2':'sw_up', # RADIAT1
                        'meantke1':'TKEs', # FLUXL: mean TKE [m2/s2] 1.5m
                        # 'h1':'shfl_s', # fluxl2: level 1 SH
                        'h2':'shfl_s', # fluxl2: level 2 SH
                        'le2':'lhfl_s', # fluxl2: level 2 LH
                        # 'meantke2':'TKE2', # FLUXL: mean TKE [m2/s2] 7.08m
                        'rawdate':'timestamp'},
                'kols': {
                        # 'taact1_avg': 'T_2M',  # RAW_TRH_PROF: TAACT1_AVG: 2m temp [°C]
                        # 'rhact1_avg': 'rh', # RAW_TRH_PROF: RHACT1_AVG: 2m rh [%]
                        'tair2':'T_2M', #	interpolated air temperature at sonic 2: 8.68m
                        'pact': 'P', #RAW: PACT: air pessure at station [hPa]
                        # 'wind_speed_4': 'VEL_10M', #SONIC_2: WIND_SPEED4: 12m wind vel [ms]
                        'meanu2':'VEL_10M', #FLUXL12: mean rotated and unfiltered u wind component (streamwise): 8.68m
                        'wind_dir2':'DIR_10M', #FLUXL12: wind direction [deg]: 8.68m
                        'meantke1':'TKEs', # FLUXL: mean TKE [m2/s2] 1.5m
                        'lw_in_avg':'lw_down', # raw
                        'lw_out_avg':'lw_up', # raw
                        'sw_in_mv_avg':'sw_down', # raw
                        'sw_out_mv_avg':'sw_up', # raw
                        # 'h1':'shfl_s', # fluxl2: level 1 SH
                        'h1':'shfl_s', # fluxl2: level 2 SH
                        # 'le1':'lhfl_s', # fluxl2: level 1 LH
                        'le1':'lhfl_s', # fluxl2: level 2 LH
                        # 'avg_wdir4': 'DIR_10M', #SONIC_2: AVG_WDIR4: 12m wind dir [deg]
                        'rawdate':'timestamp'},
                'egg': {'tair2': 'T_2M',   # RAW: Air temperature, ventilation [°C]
                        'p_avg': 'P',   # RAW: Air pressure, in loggerbox, not aerated [hPa]
                        'wind_dir2':'DIR_10M', # FLUXL: wind dir at 5.65m [deg] 
                        'meanu2':'VEL_10M', #FLUXL12
                        'h2':'shfl_s', # fluxl2: level 2 SH
                        'le2':'lhfl_s', # fluxl2: level 2 LH                           
                        'rawdate':'timestamp'},
                'weer': {'tair2': 'T_2M',   # RAW: Air temperature, ventilation [°C]
                        'p_avg': 'P',   # RAW: Air pressure, in loggerbox, not aerated [hPa]
                        'wind_dir2':'DIR_10M', # RAW: wind dir at 5.65m [deg]  
                        'meanu2':'VEL_10M', #FLUXL12                    
                        'cnr4_lw_in_wm2':'lw_down', # RADIAT1
                        'cnr4_lw_out_wm2':'lw_up', # RADIAT1
                        'cnr4_sw_in_wm2':'sw_down', # RADIAT1
                        'cnr4_sw_out_wm2':'sw_up', # RADIAT1
                        'h2':'shfl_s', # fluxl2: level 2 SH
                        'le2':'lhfl_s', # fluxl2: level 2 LH                           
                        'rawdate':'timestamp'},
                'terf': {'tair2': 'T_2M',  # RAW: Air temperature, ventilation [°C]
                        'pact': 'P', # RAW: Air pressure, in loggerbox, not aerated [hPa]
                        # 'wind_dir1':'DIR_10M', # RAW: wind dir at 6.12 m [deg]      
                        'wind_dir2':'DIR_10M', # RAW: wind dir at 11.2 m [deg]  
                        'meanu2':'VEL_10M', #FLUXL12
                        'h2':'shfl_s', # fluxl2: level 2 SH
                        'le2':'lhfl_s', # fluxl2: level 2 LH                            
                        'rawdate':'timestamp'},
                'arb': {
                        # 'taact_2m_avg': 'T_2M',  # RAW: Air temperature, ventilation [°C] 2m NOT WORKING
                        'taact_1m_avg':'T_2M',  # RAW: Air temperature, ventilation [°C] 1m
                        # 'taact_3m_avg':'T_2M',  # RAW: Air temperature, ventilation [°C] 3m
                        'pact': 'P', # RAW: Air pressure, in loggerbox, not aerated [hPa]
                        'ws_young_avg': 'VEL_10M', # RAW: Wind speed, 2m, Young Wind Monitor [ms]
                                                        # there are more for ws on1,3m but no wdir
                        'wdir_young_avg': 'DIR_10M', # RAW: Wind direction, 2m, Young Wind Monitor [deg]
                        'cnr4_lw_in_wm2':'lw_down', # RADIAT1
                        'cnr4_lw_out_wm2':'lw_up', # RADIAT1
                        'cnr4_sw_in_wm2':'sw_down', # RADIAT1
                        'cnr4_sw_out_wm2':'sw_up', # RADIAT1                          
                        'rawdate':'timestamp'},
    }

    # Read in the data
    # When downloading the data from the ACINN website, it is given as one (or
    # sometimes two) files named 'data.csv' inside a folder including the

    # station name
    name = loc.name
    # fin all files
    file_ACINN = glob2.glob(path.join(path_ACINN, f'*{name}*', 'data.csv'))
    # print number of files
    print(str(len(file_ACINN))+' files found for this station')


    if len(file_ACINN) == 0:
        raise ValueError(f'No ACINN data found for {loc.name}')

    data = pd.DataFrame()


    for file in file_ACINN: 
        data_file = pd.read_csv(file, delimiter=';', skiprows=1)
        # save data1 only if variable is insid
        for var in vars:
            # get variable nomenclature in this data/station = get key of dict value var
            var_s = getKey(acinn_v[loc.short],var.name)

            if str(var_s) in data_file.columns:
                print('found '+var.name+' as '+var_s+' in '+file)
                # extract time and var into local
                data_var=data_file[[var_s,'rawdate']]
                # converte rawdate to datetime and make it the index
                data_var.index = pd.to_datetime(data_var.rawdate)
                # cut to desired time range
                data_var = data_var[start_time: end_time]
                # add to big data
                data = pd.concat([data, data_var], ignore_index=False, axis=1)

    for c in data.columns:
        # Delete all columns that do not contain relevant data
        if c not in acinn_v[loc.short].keys():
            del(data[c])
            
        # And rename the relevant ones to be consistent
        else: #if c in acinn_v[loc.short].keys():
            data.rename(columns={c: acinn_v[loc.short][c]}, inplace=True)
    
    # drop dupliacte time columns coming from multiple source files
    data= data.T.drop_duplicates().T

    # initialize relevant variables list
    relevant_vars=[var.name for var in vars]
    relevant_vars.append("timestamp")

    # cut df to relevant variables
    data=data[relevant_vars]

    data=data.T.drop_duplicates().T
    data = data.loc[:,~data.columns.duplicated()].copy()

    print('columns: '+data.columns)
    # set_trace()

    return(data)
