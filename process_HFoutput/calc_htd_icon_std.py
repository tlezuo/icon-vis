# third party
import sys
import matplotlib.pyplot as plt
import numpy as np
import sys
from pathlib import Path
import psyplot.project as psy
import pandas as pd
import xarray
from netCDF4 import Dataset,date2num
import metpy.calc as calc
from metpy.units import units
import datetime as dt
import pandas as pd
from iconarray.plot import formatoptions # import plotting formatoptions (for use with psyplot)
import iconarray as iconvis # import self-written modules from iconarray
import io, os, sys, types
import pickle

# first party
sys.path.append('../utilities_tlezuo/')
from timefunctions import *
import varfunctions as vf
import locfunctions as lf

# define class of pc_dict
class Point_coordinates:
    def __init__(self, shortname, name, color, marker, lat, lon, altitude, iconID, iconHSURF, iconHHL, iconHFL, height_dict):
        self.short = shortname
        self.name = name
        self.color = color
        self.marker = marker
        self.lat = lat
        self.lon = lon
        self.alt = altitude
        self.iconID = iconID 
        self.iconHSURF = iconHSURF
        self.iconHHL = iconHHL
        self.iconHFL = iconHFL
        self.hdict = height_dict 


# load pc_dict and its lists
pc_dict = np.load('../utilities_tlezuo/pc_dict.npy',allow_pickle=True).item()
with open ('../utilities_tlezuo/pc_iconID_list','rb') as f:
    pc_iconID_list = pickle.load(f)
with open ('../utilities_tlezuo/pc_short_list','rb') as a:
    pc_short_list = pickle.load(a)

###############################################################################################
## DECIDE ##

# RUN
run = 'RUN3_noconv'
filepath = '/store/s83/tlezuo/'+run+'/out_std/'
savepath = '/users/tlezuo/icon-vis/data/data_std/'

# VARIABLES
# 3D
pvars_list= [vf.U, vf.V, vf.T, vf.QV,vf.P,vf.W,vf.TKE]
cvars_list = [vf.VEL, vf.DIR,vf.TD,vf.RH,vf.TH]
# 2d surface
spvars_list = [vf.T_2M,vf.QV_2M,vf.U_10M,vf.V_10M]
scvars_list = [vf.VEL_10M, vf.DIR_10M,vf.TKEs]

# TIME
startdate = dt.datetime(2019,9,13,00,00)
enddate = dt.datetime(2019,9,14,00,00)
plotfreq = '0h15min'
simdate = dt.datetime(2019,9,12,12,00) # no change,. simulation start
plotdates = pd.date_range(startdate,enddate,freq=plotfreq)

## BIG READ IN ##
# read in all nc files at once, parallelized only on node!
data=xarray.open_mfdataset(filepath+'lfff*',parallel = True)


# htd timeseries = selecting our locations
std_htd_data = data.sel(ncells=pc_iconID_list,drop=False)

# calculate new vars
# 3d
std_htd_data = std_htd_data.assign(VEL=vf.calculate_wind_vel_from_uv(std_htd_data['U'],std_htd_data['V']))
std_htd_data = std_htd_data.assign(DIR=vf.calculate_wind_dir_from_uv(std_htd_data['U'],std_htd_data['V']))
std_htd_data['P'] = std_htd_data['P']/100

# 2d
std_htd_data = std_htd_data.assign(VEL_10m=vf.calculate_wind_vel_from_uv(std_htd_data['u_10m'],std_htd_data['v_10m']))
std_htd_data = std_htd_data.assign(DIR_10m=vf.calculate_wind_dir_from_uv(std_htd_data['u_10m'],std_htd_data['v_10m']))

##  SAVE ##
std_htd_data.to_netcdf(savepath+'std_htd_data_'+run+'.nc')