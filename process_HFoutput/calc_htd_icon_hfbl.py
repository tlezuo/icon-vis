#!/usr/bin/env python3

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

###############################################################################################
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
filepath = '/store/s83/tlezuo/'+run+'/out_hfbl/'
savepath = '/users/tlezuo/icon-vis/data/data_hfbl/'

###############################################################################################
# VARIABLES
pvars_list= [vf.T, vf.TKEVELtend, vf.Ttend_clcov, vf.Ttend_drag, vf.Ttend_pconv, vf.Ttend_radlw, vf.Ttend_radsw, vf.Ttend_turb, vf.Ttend_dyn] # 3d

# TIME
startdate = dt.datetime(2019,9,13,00,00)
enddate = dt.datetime(2019,9,14,00,00)
plotfreq = '0h0min10s'
simdate = dt.datetime(2019,9,12,12,00) # no change,. simulation start
plotdates = pd.date_range(startdate,enddate,freq=plotfreq)

###############################################################################################
## BIG READ IN ##

# read in all nc files at once, parallelized only on node!
data=xarray.open_mfdataset(filepath+'lfffhfbl*', parallel=True)

# htd timeseries = selecting our locations
hfbl_htd_data = data.sel(ncells=pc_iconID_list,drop=False)

# surface timeseries
hfbl_ts_data = hfbl_htd_data.sel(height=80, height_2 = 80)

# surface timeseries integrated: sum of tendencies each hour *10s 
sumbl_ts_data=hfbl_ts_data.groupby("time.hour").sum(dim='time')*10

# htd timeseries integrated
# only for 3 rs and lidar locations
pc_short_list_hfblvert = ['ifl','kols','murs']
pc_iconID_list_hfblvert = [17,12,18] # these are the new indices in the subset (pc_dict[kols].subsetID = 12)
hfbl_htd_data_hfblvert = hfbl_htd_data.sel(ncells=pc_iconID_list_hfblvert,drop=False)
sumbl_htd_data=hfbl_htd_data_hfblvert.groupby("time.hour").sum(dim='time')*10

###############################################################################################
## SAVE ##
# hfbl_htd_data.to_netcdf(savepath+'hfbl_htd_data_'+run+'.nc')
# hfbl_ts_data.to_netcdf(savepath+'hfbl_ts_data_'+run+'.nc')
sumbl_ts_data.to_netcdf(savepath+'sumbl_htd_data_'+run+'.nc')
sumbl_htd_data.to_netcdf(savepath+'sumbl_htd_data_'+run+'.nc')