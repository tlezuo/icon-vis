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

# 2) vertical cross sections
class VC_coordinates:
    def __init__(self, name, short, dist, latmin, lonmin, namemin, latmax, lonmax, namemax):
        self.name = name
        self.short = short
        self.dist = dist
        self.latmin = latmin
        self.lonmin = lonmin
        self.namemin = namemin
        self.latmax = latmax
        self.lonmax = lonmax
        self.namemax = namemax
    

####################################################################
## COORDINATE SET ##
# 1) local coordinates for points
# pc_dict = np.load('pc_dict.npy',allow_pickle='TRUE').item()

# 2) area coordinates for hor cross
# 2)a) icon domain
icon_domain = Area_coordinates()
icon_domain.name = 'icon_domain'
icon_domain.latmin = 42.7
icon_domain.lonmin = 1
icon_domain.latmax = 49.7
icon_domain.lonmax = 16.3
icon_domain.locmarks = []
icon_domain.locmarks.extend(['inn','mun','ven'])


# 2)b) eastern alps
eastern_alps = Area_coordinates()
eastern_alps.name = 'eastern_alps'
eastern_alps.latmin = 45
eastern_alps.lonmin = 9
eastern_alps.latmax = 49
eastern_alps.lonmax = 15
eastern_alps.locmarks = []
eastern_alps.locmarks.extend(['inn','mun','ven','sal','bz','ver'])

# 2)c) innsbruck area
inn_area = Area_coordinates()
inn_area.name = 'inn_area'
inn_area.latmin = 46.75
inn_area.lonmin = 10
inn_area.latmax = 48.3
inn_area.lonmax = 12.7
inn_area.locmarks = []
inn_area.locmarks.extend(['inn', 'hoch', 'kols', 'weer', 'egg', 'weer', 'terf', 'arb', 'iun', 'ifl'])

# 2)c) local area
local_area = Area_coordinates()
local_area.name = 'local_area'
local_area.latmin = 47.1
local_area.lonmin = 11
local_area.latmax = 47.5
local_area.lonmax = 12
local_area.locmarks = []
local_area.locmarks.extend(['hoch', 'kols', 'weer', 'egg', 'weer', 'terf', 'arb', 'iun', 'ifl'])

# 2)c) domain visualisation
domain_vis = Area_coordinates()
domain_vis.name = 'domain_vis'
domain_vis.latmin = 40
domain_vis.lonmin = 0
domain_vis.latmax = 50
domain_vis.lonmax = 20
domain_vis.locmarks = []
domain_vis.locmarks.extend(['inn'])

# 2)c) super zoom ifl
zoom = Area_coordinates()
zoom.name = 'zoom'
zoom.latmin = 47.2
zoom.lonmin = 11.2
zoom.latmax = 47.4
zoom.lonmax = 11.5
zoom.locmarks = []
zoom.locmarks.extend(['ifl'])

# 2)c) local local area: only stations visible
stat_only = Area_coordinates()
stat_only.name = 'stations area'
stat_only.latmin = 47.2
stat_only.lonmin = 11.22
stat_only.latmax = 47.38
stat_only.lonmax = 11.81
stat_only.locmarks = []
stat_only.locmarks.extend(['hoch', 'kols', 'weer', 'egg', 'weer', 'terf', 'arb', 'iun', 'ifl'])

# 2)c) local local area: only stations visible
nies = Area_coordinates()
nies.name = 'niesen'
nies.latmin = 46.59
nies.lonmin = 7.58
nies.latmax = 46.69
nies.lonmax = 7.71
nies.locmarks = []

# 2) vertical cross sections
# MY OWN
VCS_kols = VC_coordinates('VCS at Kolsass' , 'VCS_KOL', 45,
                         47.50, 11.50, 'Simmering',
                         47.095102, 11.721242, 'S', )
VCS_alp = VC_coordinates('VCS at Alpbach' , 'VCS_ALP', 45,
                         47.643407, 11.743848, 'Kreuth',
                         47.219784, 12.180235,'Krimml', )
VCS_hai = VC_coordinates('VCS at Haiming' , 'VCS_HAI', 45,
                         47.418938, 10.989080, 'Zugspitze',
                         47.034280, 11.134006,'Stubaier Gletscher', )
VCS_kuf = VC_coordinates('VCS at Kufstein' , 'VCS_KUF', 45,
                         47.657884, 11.944264, 'Ruchenkopf',
                         47.516713, 12.426752,'St Johann', )
VCS_border = VC_coordinates('VCS at border' , 'VCS_Border', 45,
                         47.720323, 11.867407, 'Schielersee',
                         47.730065, 12.463968,'Unterwössen', )

# from Toni
VCS_up_valley = VC_coordinates('VCS in upper Inn Valley' , 'VCS_up', 100,
                         47.46831375, 11.46113725, '__',
                         47.09595825,11.52701275,'__', )
VCS_down_valley = VC_coordinates('VCS in lower Inn Valley' , 'VCS_down', 100,
                         47.499804, 11.642547, '__',
                         47.227952,11.802475,'__', )

## LIDARS VCS, from Toni
VCS_lidars1 = VC_coordinates('Cross Valley Section along Lidar VCS' , 'VCS_lidars1', 100,
                         47.33668, 11.6024, 'Eggen',
                         47.305290,11.622231,'Weer', )
VCS_lidars2 = VC_coordinates('Cross Valley Section along Lidar VCS' , 'VCS_lidars2', 100,
                         47.305290,11.622231,'Weer', 
                         47.27108, 11.63841, 'Hoch' )
VCS_lidars_full = VC_coordinates('VCS Lidar' , 'VCS_lidars_full', 100,
                          47.33668, 11.6024, 'Eggen', 
                         47.27108, 11.63841, 'Hoch' )

## ALONG VALLEY, from TONI
VCS_av1 = VC_coordinates('VCS along valley' , 'VCS_av1', 100,
                        47.311373, 11.063994, 'A', 
                        47.241373,  11.383994, 'B' )
VCS_av2 = VC_coordinates('VCS along valley' , 'VCS_av2', 100,
                        47.241373,  11.383994, 'B',
                        47.305290, 11.622231, 'Kols', )
VCS_av3 = VC_coordinates('VCS along valley' , 'VCS_av3', 100,
                        47.305290, 11.622231, 'Kols', 
                        47.425023,  11.843750, 'C',)
VCS_av4 = VC_coordinates('VCS along valley' , 'VCS_av4', 100,
                        47.425023,  11.843750, 'C',
                        47.515023,  12.083750, 'D', )
VCS_av5 = VC_coordinates('VCS along valley' , 'VCS_av5', 100,
                        47.515023,  12.083750, 'D', 
                        47.635023,  12.193750, 'E',)
VCS_av6 = VC_coordinates('VCS along valley' , 'VCS_av6', 100,
                        47.635023,  12.193750, 'E',
                        47.85023, 12.103750, 'F', )
VCS_av_full = VC_coordinates('VCS along valley' , 'VCS_av_full', 100,
                        47.311373, 11.063994, 'A', 
                        47.85023, 12.103750, 'F', )
VCS_av_full_T = VC_coordinates('VCS along valley' , 'VCS_AV', 100,
                        47.241373,  11.383994, 'B', 
                        47.515023,  12.083750, 'D', )

## ALPINE PUMPING, from Toni
VCS_ap1 = VC_coordinates('VCS N-S alpine pumping' , 'VCS_ap1', 100,
                        47.095102, 11.721242, 'S', 
                        47.305290, 11.622231, 'Kols', )
VCS_ap2 = VC_coordinates('VCS N-S alpine pumping' , 'VCS_ap2', 100,
                        47.305290, 11.622231, 'Kols', 
                        47.861014, 11.271361, 'N', )
VCS_ap_full = VC_coordinates('VCS N-S alpine pumping' , 'VCS_ap_full', 100,
                        47.095102, 11.721242, 'S', 
                        47.861014, 11.271361, 'N', )

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


