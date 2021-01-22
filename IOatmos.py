import time
from datetime import datetime, timedelta
import os, sys, string
from netCDF4 import Dataset
import numpy as np

"""
Created by Trond Kristiansen
https://github.com/trondkr/model2roms
"""

def help ():
    """
    This function generates the initial netcdf atmospheric forcing file for the U and V wind component
    for ROMS. Methods:

    def createNetCDFFileUV(grdROMS, outfilename, myformat):
    """

def createNetCDFFileUV(grdROMS, outfilename, myformat, mytype):
    if (myformat=='NETCDF4'):
        myzlib = True
    else:
        myzlib = False

    if os.path.exists(outfilename):
        os.remove(outfilename)

    f1             = Dataset(outfilename, mode='w', format=myformat)
    f1.title       = "Atmospheric forcing file"
    f1.grdFile     = "%s"%(grdROMS.grdfilename)
    f1.history     = 'Created ' + time.ctime(time.time())
    f1.source      =  "{} ({})".format(confM2R.authorname, confM2R.authoremail)
    f1.Conventions = "CF-1.0"

    """ Define dimensions """
    f1.createDimension('xi_rho',  grdROMS.xi_rho)
    f1.createDimension('eta_rho', grdROMS.eta_rho)
    f1.createDimension('wind_time', None)
    
    vnc               = f1.createVariable('lon_rho', 'd', ('eta_rho','xi_rho',),zlib=myzlib, fill_value=grdROMS.fill_value)
    vnc.long_name     = 'Longitude of RHO-points'
    vnc.units         = 'degree_east'
    vnc.standard_name = 'longitude'
    vnc[:,:]          = grdROMS.lon_rho

    vnc               = f1.createVariable('lat_rho', 'd', ('eta_rho','xi_rho',),zlib=myzlib, fill_value=grdROMS.fill_value)
    vnc.long_name     = 'Latitude of RHO-points'
    vnc.units         = 'degree_north'
    vnc.standard_name = 'latitude'
    vnc[:,:]          = grdROMS.lat_rho

    v_time            = f1.createVariable('wind_time', 'd', ('wind_time',),zlib=myzlib, fill_value=grdROMS.fill_value)
    v_time.long_name  = 'Days since 1948-01-01 00:00:00'
    v_time.units      = 'Days since 1948-01-01 00:00:00'
    v_time.field      = 'time, scalar, series'
    v_time.calendar   = 'standard'

    v_temp_west               = f1.createVariable('Vwind', 'f', ('wind_time', 'eta_rho', 'xi_rho',),zlib=myzlib, fill_value=grdROMS.fill_value)
    v_temp_west.long_name     = "Eta-component of wind"
    v_temp_west.units         = "meter second-1"
    v_temp_west.field         = "Vwind, scalar, series"
    v_temp_west.missing_value = grdROMS.fill_value
    v_temp_west.time          = "wind_time"

    v_temp_west               = f1.createVariable('Uwind', 'f', ('wind_time', 'eta_rho', 'xi_rho',),zlib=myzlib, fill_value=grdROMS.fill_value)
    v_temp_west.long_name     = "Xi-component of wind"
    v_temp_west.units         = "meter second-1"
    v_temp_west.field         = "Uwind, scalar, series"
    v_temp_west.missing_value = grdROMS.fill_value
    v_temp_west.time          = "wind_time"

    f1.close()
