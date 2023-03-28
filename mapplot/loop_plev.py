# Load modules
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import cmcrameri.cm as cmc
import cartopy.feature as cf                                                                                                        
from pathlib import Path
import psyplot.project as psy
import sys
import datetime as dt
from iconarray.plot import formatoptions # import plotting formatoptions (for use with psyplot)
import iconarray as iconvis # import self-written modules from iconarray
import matplotlib.patheffects as PathEffects

# #own scripts
sys.path.append('../utilities_tlezuo/')
import timefunctions as tf
import varfunctions as vf
import locfunctions as lf

#  define class of pc_dict
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

###################################### GLOBAL #######################################
simdate = dt.datetime(2019,9,12,12,00)
grid_file = '../data/example_data/grids/ICON-1E_DOM01.nc'

###################################### DECIDE #######################################
area=lf.lower_valley
pc_innvalley_list = ['ALP','RIN','UNI','kols']
# plevs_list = [70000., 75000., 80000., 85000., 90000.,92500., 95000.]
plevs_list = [85000., 95000.]
startdate_model = dt.datetime(2019,9,12,12,00)
enddate_model = dt.datetime(2019,9,14,12,00)
pdates_list = pd.date_range(startdate_model,enddate_model,freq='1H')
lake_tf=True
bord_tf=True
ab='a)'
x_grid= False
y_grid = False
plottype = 'poly' # 'poly' = triangles,'contourf'= smoothed, 'contour'= lines
# variable
pvar = vf.T


###################################### LOOP 1 TIME ####################################
for pdate in pdates_list:
    lt = tf.get_lt(pdate,simdate)
    # print(lt)
    lfff_name =tf.lfff_name(lt)
    filename = 'lfffp'+lfff_name[4:]
    nc_file = '/store/s83/tlezuo/RUN5_extended/out_p/'+filename
    data = psy.open_dataset(nc_file)
    print(filename)
    
###################################### LOOP 2 plev ####################################
    for PL in plevs_list: 
        
        # get cmap levels
        if PL == 85000:
            levels = np.arange(282,288,0.5)
        elif PL == 95000:
            levels = np.arange(287,295,0.5)

        # path sepcific to area and pressure level
        plotpath = 'plots'+str(area.name)+'/'+str(PL)+'/'
        Path(plotpath).mkdir(parents=True, exist_ok=True)
        # title specific to area, variable, time
        title = 'temperature and wind at '+pdate.strftime('%d.%m at %H%M UTC')+' on '+str(PL/100)+' hPa'
        plotname = str(area.name)+'_'+pvar.name+'_'+str(PL)+'_'+pdate.strftime('%y%m%d_%H%M')+'.png'

        ###############################################################################################
        ## PLOT ##
        # psy.plot.mapcombined(data.sel(plev=PL),
        #                                 name=[['T', 
        #                                 ['U', 'V']]], 
        #                                 map_extent = [area.lonmin, area.lonmax, area.latmin, area.latmax],
        #                                 plot=plottype,
        #                                 arrowsize=200,
        #                                 bounds = np.arange(285,295,0.5),
        #                                 title = title,
        #                                 cmap=pvar.cmap,
        #                                 cticks=np.arange(285,295,0.5),
        #                                 clabel = pvar.name+' '+pvar.units,
        #                                 xgrid = x_grid, ygrid = y_grid,
        #                                 cbar = 'b',#'b,r False',
        #                                 lakes=lake_tf,
        #                                 borders=bord_tf,
        #                                 # background = 'grey'
        #                                 )
        psy.plot.mapplot(data.sel(plev=PL),
                                name=pvar.name, 
                                map_extent = [area.lonmin, area.lonmax, area.latmin, area.latmax],
                                plot=plottype,
                                bounds =levels,
                                title = title,
                                cmap='coolwarm',
                                cticks=levels,
                                clabel = pvar.name+' '+pvar.units,
                                xgrid = x_grid, ygrid = y_grid,
                                cbar = 'b',#'b,r False',
                                lakes=lake_tf,
                                borders=bord_tf,
                                # background = 'grey'
                                  
                                  )

        # figure
        fig = plt.figure
        fig = plt.gcf()
        fig.set_size_inches(10,10)
        plt.rcParams['figure.dpi'] = 300
        # plt.title(ab,loc='left')

        #annotations
        pc_dict['murs'].short = 'MUA'
        pc_dict['kols'].short = 'KOL'
        pc_dict['egg'].short = 'EGG'
        pc_dict['terf'].short = 'TER'
        pc_dict['arb'].short = 'ARB'
        pc_dict['hoch'].short = 'HOC'
        pc_dict['weer'].short = 'WEE'
        for locmark in pc_innvalley_list:
            
            if pc_dict[locmark].short == 'KOL' or pc_dict[locmark].short == 'MUA' :
                col_marker = 'w'#'violet'
            else:
                col_marker = 'k'#'violet'
            pos_lon, pos_lat = iconvis.add_coordinates(pc_dict[locmark].lon,pc_dict[locmark].lat,area.lonmin,area.lonmax,area.latmin,area.latmax)
            mark = fig.axes[0].plot(pos_lon, pos_lat,color=col_marker,marker='.',markersize=5, transform=fig.axes[0].transAxes,label=pc_dict[locmark].name,linestyle='None',zorder=101)
                                    #,path_effects=[PathEffects.Stroke(linewidth=3, foreground='k'), PathEffects.Normal()]) 
            mtxt = fig.axes[0].annotate(text=pc_dict[locmark].short,xy=(pos_lon, pos_lat),xycoords='axes fraction',zorder=100,color=col_marker, style = 'oblique',weight='bold')

        # save figure
        plt.savefig(plotpath+plotname,dpi=300)