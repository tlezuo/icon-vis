import sys
from pathlib import Path
import numpy as np
import xarray as xr


####################################################################
## CLASS DEFINITION ##
# 1) area coordinates for hor cross
class Area_coordinates:
    # need 2 locations to delimit area
    name = 'AREANAME'
    latmin = 0
    lonmin = 0
    latmax = 0
    lonmax = 0
    locmarks = []

# 2) local coordinates for points
class Point_coordinates:
    def __init__(self, shortname, name, color, marker, lat, lon, altitude, iconID):
        self.short = shortname
        self.name = name
        self.color = color
        self.marker = marker
        self.lat = lat
        self.lon = lon
        self.alt = altitude
        self.iconID = iconID # gridpoint used in icon for this location

# 3) station metadata
class Station_metadata:
    def __init__(self, shortname, name, lat, lon, altitude, altitude_icon, height_dict):
        self.short = shortname
        self.name = name
        self.lat = lat
        self.lon = lon
        self.alt = altitude
        self.iconalt = altitude_icon
        self.hdict = height_dict # gridpoint used in icon for this location

####################################################################
## COORDINATE SET ##
# 1) local coordinates for points
# 1) a) Innsbruck
inn = Point_coordinates('inn','Innsbruck', 'k','.',47.26,11.42,0,0)
# Salzburg
sal = Point_coordinates('sal','Salzburg', 'k','.', 47.79,13.09,0,0)
# München
mun = Point_coordinates('mun','Munich','k','.',48.14,11.57,0,0)
# Bozen
bz = Point_coordinates('bz','Bozen','k','.', 46.49,11.36,0,0)
# Verona
ver = Point_coordinates('ver','Verona', 'k','.',45.43,10.99,0,0)
# Venice
ven = Point_coordinates('ven','Venice','k','.',45.33,12.27,0,0)

## meteo stations and rs
# iun
iun = Point_coordinates('iun','stat Universität','gold',"$u$",47.2642889,11.3861614,0,0)
# ifl station + rs
ifl = Point_coordinates('ifl','Station & Radiosounding Innsbruck Airport','gold',"$f$",47.25846,11.3521825,579,808260)
ifl_stat = Station_metadata('ifl','Station & Radiosounding Innsbruck Airport',47.25846,11.3521825,579,206,{})
# mu rs
murs = Point_coordinates('murs','Radiosounding Munich','gold',".",48.25,11.55,492,753969)
murs_stat = Station_metadata('murs','Radiosounding Munich',48.25,11.55,492,206,{})

## ibox stations
# Hoch
hoch = Point_coordinates('hoch','Hochhaueser','violet',"$h$",47.28755, 11.63122,1009,810547)
hoch_stat = Station_metadata('hoch','Hochhaueser',47.28755, 11.63122,1009,1114,{'T_2M':6.83,
                                                                            'VEL_10M':7.08,
                                                                            'DIR_10M':7.08,
                                                                            'TKEs':1.50,
                                                                            'turb2':7.08,
                                                                            'rad':2.00,
                                                                            })

# Kols
kols = Point_coordinates('kols','Kolsass','darkviolet',"$k$",47.305, 11.622,547,810552)
kols_stat = Station_metadata('kols','Kolsass',47.305, 11.622,547,551,{'T_2M':8.68,#2,4,8.68,16.93
                                                                            'VEL_10M':8.68, #2,4,6,12
                                                                            'DIR_10M':8.68,
                                                                            'TKEs':4.00,
                                                                            'turb2':8.68,
                                                                            'turb3':16.93,
                                                                            'rad':2.00,
                                                                            })

# Eggen
egg = Point_coordinates('egg','Eggen','violet',"$e$",47.3165,11.6162,892,810583)
# Weer
weer = Point_coordinates('weer','Weerberg','violet',"$w$",47.299,11.672,930,810560)
# Terf
terf = Point_coordinates('terf','Terfens','violet',"$t$", 47.325538,11.65247,575,810571)
# Arb
arb = Point_coordinates('arb','Arbeser','violet',"$a$",47.320654,11.746592,2015,811094)


# 2) area coordinates for hor cross
# 2)a) icon domain
icon_domain = Area_coordinates()
icon_domain.name = 'icon_domain'
icon_domain.latmin = 42.7
icon_domain.lonmin = 1
icon_domain.latmax = 49.7
icon_domain.lonmax = 16.3
icon_domain.locmarks = []
icon_domain.locmarks.extend([inn,mun,ven])


# 2)b) eastern alps
eastern_alps = Area_coordinates()
eastern_alps.name = 'eastern_alps'
eastern_alps.latmin = 45
eastern_alps.lonmin = 9
eastern_alps.latmax = 49
eastern_alps.lonmax = 15
eastern_alps.locmarks = []
eastern_alps.locmarks.extend([inn,mun,ven,sal,bz,ver])

# 2)c) innsbruck area
inn_area = Area_coordinates()
inn_area.name = 'inn_area'
inn_area.latmin = 46.5
inn_area.lonmin = 10
inn_area.latmax = 48
inn_area.lonmax = 13
inn_area.locmarks = []
inn_area.locmarks.extend([inn, hoch, kols, weer, egg, weer, terf, arb, iun, ifl])

# 2)c) local area
local_area = Area_coordinates()
local_area.name = 'local_area'
local_area.latmin = 47.1
local_area.lonmin = 11
local_area.latmax = 47.5
local_area.lonmax = 12
local_area.locmarks = []
local_area.locmarks.extend([hoch, kols, weer, egg, weer, terf, arb, iun, ifl])

# 2)c) domain visualisation
domain_vis = Area_coordinates()
domain_vis.name = 'domain_vis'
domain_vis.latmin = 40
domain_vis.lonmin = 0
domain_vis.latmax = 50
domain_vis.lonmax = 20
domain_vis.locmarks = []
domain_vis.locmarks.extend([inn])

# 2)c) super zoom ifl
zoom = Area_coordinates()
zoom.name = 'zoom'
zoom.latmin = 47.2
zoom.lonmin = 11.2
zoom.latmax = 47.4
zoom.lonmax = 11.5
zoom.locmarks = []
zoom.locmarks.extend([ifl])

# 2)c) local local area: only stations visible
stat_only = Area_coordinates()
stat_only.name = 'stations area'
stat_only.latmin = 47.2
stat_only.lonmin = 11.22
stat_only.latmax = 47.38
stat_only.lonmax = 11.81
stat_only.locmarks = []
stat_only.locmarks.extend([hoch, kols, weer, egg, weer, terf, arb, iun, ifl])

# 2)c) local local area: only stations visible
nies = Area_coordinates()
nies.name = 'niesen'
nies.latmin = 46.59
nies.lonmin = 7.58
nies.latmax = 46.69
nies.lonmax = 7.71
nies.locmarks = []




####################################################################
## COORDINATE SET ##
# print(eastern_alps.locmarks)
# for mark in eastern_alps.locmarks:
#     print(mark.name)

def ind_from_latlon(lats, lons, lat, lon):
    """Find the nearest neighbouring index to given location.
    Args:
        lats (2d array):            Latitude grid
        lons (2d array):            Longitude grid
        lat (float):                Latitude of location
        lon (float):                Longitude of location
    Returns:
        int     Index of nearest grid point.
    """
    dist = [
        np.sqrt((lats[i] - lat) ** 2 + (lons[i] - lon) ** 2) for i in range(len(lats))
    ]
    ind = np.where(dist == np.min(dist))[0][0]

    return ind

def get_grid_names(ds_grid, verbose=False):

    possible_lats_names = ["clat", "clat_1"]
    possible_lons_names = ["clon", "clon_1"]
    possible_height_names = [
        "HHL",
        "hhl",
        "HEIGHT",
    ]
    possible_index_names = ["cells_1", "ncells", "cells"]

    for coo in ds_grid.coords:

        # clat
        for name in possible_lats_names:
            if name == coo:
                lats_name = name
                break

        # clon
        for name in possible_lons_names:
            if name == coo:
                lons_name = name
                break

    for var in ds_grid.variables:

        # half level height
        for name in possible_height_names:
            if name == var:
                height_name = name
                break

    for dim in ds_grid[height_name].dims:

        # cell index name for height variable
        for name in possible_index_names:
            if name == dim:
                height_index_name = name
                break

    return lats_name, lons_name, height_name, height_index_name

def index_height_from_height_file(lat, lon, grid):
    """Retrieve index and height for specific grid point.
    Args:
        lat (float): latitude
        lon (float): longitude
        grid (str): grid file (netcdf)
    Returns:
        index (int)
        height (1-dimensional np array)
        size (int): size of grid file
    """

    # load grid file
    if Path(grid).is_file():
        ds_grid = xr.open_dataset(grid).squeeze()
    else:
        print("Grid file does not exist!")
        sys.exit(1)

    lats_name, lons_name, height_name, height_index_name = get_grid_names(ds_grid)
    # load latitude and longitude grid of constants file
    lats_grid = ds_grid[lats_name].values
    lons_grid = ds_grid[lons_name].values

    # convert from radians to degrees if given in radians
    if lats_grid.max() < 2.0 and lons_grid.max() < 2.0:
        lats_grid = np.rad2deg(lats_grid)
        lons_grid = np.rad2deg(lons_grid)

    # find index closest to specified lat, lon (in grid file)
    ind = ind_from_latlon(lats_grid, lons_grid, lat, lon, False)

    # load HEIGHT from grid file
    ds_height = ds_grid[height_name]
    height = ds_height.isel(**{height_index_name: ind}).values

    return ind, height, lats_grid.size


