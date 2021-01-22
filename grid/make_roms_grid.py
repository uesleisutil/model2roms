#!/home/uesleisutil/anaconda3/bin/python python
# -*- coding: utf-8 -*-

"""
File name:     make_roms_grid.py
Author:        Ueslei Adriano Sutil
Email:         uesleisutil1@gmail.com
Created:       17 June 2016
Last modified: 04 July 2019
Version:       6.1

Variable options:
    - grd_name         As mentioned in gridid.txt Pyroms file;
    - intrp_method    'linear' or 'nn' for nearest neightborhood;
    - grid_resolution  The ETOPO1 spatial resolution is 1/60°. If you desire a 1/12° spatial resolution, choose 5 because 5/60 = 1/12°
                       The SRTM30+ spatial resolution is 1/120°. If you desire a 1/12° spatial resolution, choose 10 because 10/120 = 1/12°;
    - max_depth        The grid max depth (m). 

Download ETOPO1 here:
    https://rda.ucar.edu/datasets/ds759.4/index.html#!description
Reference: 
    Amante, C. and B. W. Eakins, 2009: ETOPO1 1 Arc-Minute Global Relief Model: Procedures, Data Sources and Analysis. 24, NOAA Technical Memorandum NESDIS NGDC, 19 pp.

Download SRTM30_plus here:
    https://topex.ucsd.edu/WWW_html/srtm30_plus.html
Reference:
    Becker, J. J., D. T. Sandwell, W. H. F. Smith, J. Braud, B. Binder, J. Depner, D. Fabre, J. Factor, S. Ingalls, S-H. Kim, R. Ladner, K. Marks, S. Nelson, A. Pharaoh, R. Trimmer, J. Von Rosenberg, G. Wallace, P. Weatherall., Global Bathymetry and Elevation Data at 30 Arc Seconds Resolution: SRTM30_PLUS, Marine Geodesy, 32:4, 355-371, 2009.
"""

from   mpl_toolkits.basemap import Basemap
from   bathy_smoother       import bathy_tools, bathy_smoothing
import mpl_toolkits.basemap as mp
import numpy                as np
import netCDF4
import pyroms
import pyroms_toolbox
from   sty                 import bg

# Outputs names.
grd_name              = 'atlsw_op'
grd_final             = grd_name+'_grd.nc'
etopo1_dir            = '/home/ueslei/Documents/model2roms/grid/etopo1.nc'
srtm_dir              = '/home/ueslei/Documents/model2roms/grid/topo30_atlsw.nc'

# Grid settings.
hmin                  = 20
theta_b               = 0.6
theta_s               = 5.0
Tcline                = 50
N                     = 30
rmax                  = 0.45
intrp_method          = 'linear'
grid_resolution       = 10
max_depth             = -5000


# NOTE: Do not need to change above.

# Open bathymetry file.
print(bg.da_cyan+"Choose the Digital Elevation Model: (1) ETOPO1 or (2) SRTM30+, then press the [ENTER] buttom:"+bg.rs)
grid_choice = input()
if grid_choice=="1":
    data = netCDF4.Dataset(etopo1_dir, 'r')
    lons = data.variables[u'x'][:]
    lats = data.variables[u'y'][:]
    cota = data.variables[u'z'][:]
    cota = np.array(cota, dtype='float32')
elif grid_choice=="2":
    data = netCDF4.Dataset(srtm_dir, 'r')
    lons = data.variables[u'lon'][:]
    lats = data.variables[u'lat'][:]
    cota = data.variables[u'z'][:]
    cota = np.array(cota, dtype='float32')

# Grid max depth.
hcopy = cota.copy()
for i in range(len(cota[:,1])):
    for j in range(len(cota[1,:])):
        if hcopy[i,j]<=max_depth:
            hcopy[i,j]=max_depth
cota = hcopy

# Generating the desired resolution.
resol     = lons[1]-lons[0]
resol     = (resol*grid_resolution)
print(bg.da_cyan+'The ROMS grid spatial resolution is:',resol*100,'km x',resol*100,'km.'+bg.rs)
lons1     = np.arange(lons.min(),lons.max(),resol)
lats1     = np.arange(lats.min(),lats.max(),resol)
lon1,lat1 = np.meshgrid(lons1, lats1)
cota1     = mp.interp(cota,lons,lats,lon1,lat1,checkbounds=False,masked=False,order=1)

# Grid dimension.
Mm   = len(lats1)
Lm   = len(lons1)
lon0 = lons1.min() ; lat0 =lats1.max()
lon1 = lons1.min() ; lat1 =lats1.min()
lon2 = lons1.max() ; lat2 =lats1.min()
lon3 = lons1.max() ; lat3 =lats1.max()

# Choosing grid projection.
print(bg.da_cyan+"Choose the ROMS grid projection: (1) Mercator, (2) Polar or (3) Lambert Conformal, then press the [ENTER] buttom:"+bg.rs)
map_projection = input()
if map_projection=='1':
    map = Basemap(projection='merc', llcrnrlon=lons1.min(), llcrnrlat=lats1.min(),urcrnrlon=lons1.max(), urcrnrlat=lats1.max(), resolution='f')
if map_projection=='2':
    map = Basemap(projection='spstere',boundinglat=-45,lon_0=90,resolution='f') 
if map_projection=='3':
    map = Basemap(width=12000000,height=9000000,rsphere=(6378137.00,6356752.3142),resolution='f',area_thresh=1000.,projection='lcc',lat_0=lat0, lat_1=lat1, lat_2=lat3, lon_0 =lon0)

lonp = np.array([lon0, lon1, lon2, lon3])
latp = np.array([lat0, lat1, lat2, lat3])
beta = np.array([1., 1., 1., 1.])

# Start generating the new grid.
hgrd       = pyroms.grid.Gridgen(lonp, latp, beta, (Mm,Lm),proj=map)
lonv, latv = map(hgrd.x_vert, hgrd.y_vert, inverse=True)
hgrd       = pyroms.grid.CGrid_geo(lonv, latv, map)

for verts in map.coastsegs:
    hgrd.mask_polygon(verts)

# Check the ocean-continent mask and edit if necessary.
coast = pyroms.utility.get_coast_from_map(map)
pyroms.grid.edit_mask_mesh_ij(hgrd, coast=coast)

# Interpolate new bathymetry.
h = mp.interp(cota1,lons1,lats1,hgrd.lon_rho,hgrd.lat_rho,checkbounds=False,masked=False,order=0)

# Save raw bathymetry.
hraw=h.copy()

# ROMS depth is positive.
hh = -h

# Depth deeper than hmin.
h = np.where(hh < hmin, hmin, hh)

# Smooth the raw bathy using the direct iterative method from Martinho and Batteen (2006).
RoughMat = bathy_tools.RoughnessMatrix(h, hgrd.mask_rho)
print ('Currently, the max roughness value is: ', RoughMat.max())
h        = bathy_smoothing.smoothing_Positive_rx0(hgrd.mask_rho, h, rmax)
h        = bathy_smoothing.smoothing_Laplacian_rx0(hgrd.mask_rho, h, rmax)
RoughMat = bathy_tools.RoughnessMatrix(h, hgrd.mask_rho)
print (('After the filters, the max roughness value is: ', RoughMat.max()))
hgrd.h   = h

# Vertical levels.
vgrd = pyroms.vgrid.s_coordinate_4(h, theta_b, theta_s, Tcline, N, hraw=hraw)

# ROMS grid.
grd = pyroms.grid.ROMS_Grid(grd_name, hgrd, vgrd)

# Write grid to netcdf file.
pyroms.grid.write_ROMS_grid(grd, grd_final)
