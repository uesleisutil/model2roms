import time, calendar
from netCDF4 import Dataset, date2num, num2date
import model2roms
import clim2bry
import grd
import numpy as np
import atmosForcing
import sys
from datetime import datetime, timedelta
import os

"""
Created by Trond Kristiansen
https://github.com/trondkr/model2roms
"""

class Model2romsConfig(object):

    # Define abbreviation for the run: sued to name output files etc.
    def defineabbreviation(self):
        return {"antarctic": "antarctic"}[self.outgrid]

    def showinfo(self):
        print('Creating condition files for ROMS from %s/%s to %s/%s' % (self.start_year,self.start_month, self.end_year,self.end_month))

    def formatdatesforoutputnames(self):
        # Format the date for use in output filenames
        startmonth = ("0%s" % self.start_month if self.start_month < 10 else "%s" % self.start_month)
        startday = ("0%s" % self.start_day if self.start_day < 10 else "%s" % self.start_day)
        endmonth = ("0%s" % self.end_month if self.end_month < 10 else "%s" % self.end_month)
        endday = ("0%s" % self.end_day if self.end_day < 10 else "%s" % self.end_day)

        modelperiod = str(self.start_year) + str(startmonth) + str(startday) + '_to_' + str(self.end_year) + str(
            endmonth) + str(
            endday)

        return modelperiod

    def defineoutputfilenames(self):
        # Get string representation of start and end dates
        modelperiod = self.formatdatesforoutputnames()

        # Name of output files for CLIM, BRY, and INIT files
        climname = self.abbreviation + '_cli.nc'
        initname = self.abbreviation + '_ini.nc'
        bryname = self.abbreviation + '_bry.nc'

        return climname, initname, bryname

    # Define the global variables to be used for each type of input data. Not all input datasets contains information on
    # e.g. sea ice so those variables can not be included. the SODA3si for example does not contain ssh.
    # OPTIONS: ['temperature', 'salinity', 'ssh', 'uvel', 'vvel', 'ageice', 'uice', 'vice', 'aice', 'hice']

    def defineglobalvarnames(self):

        return {'SODA3': ['temperature', 'salinity', 'ssh', 'uvel', 'vvel','uice','vice','aice','hice','snow_thick','ageice'],
                'GLORYS': ['temperature', 'salinity', 'ssh', 'uvel', 'vvel', 'uice', 'vice', 'aice', 'snow_thick']}[self.indatatype]
        #return {'SODA3': ['temperature', 'salinity', 'ssh', 'uvel', 'vvel','uice','vice','aice','hice','snow_thick','ageice'],
        #        'GLORYS': ['temperature', 'salinity', 'ssh', 'uvel', 'vvel']}[self.indatatype]

    # Define the corresponding name of the variables in the input dataset files. This list needs to correspond
    # exactly with the list given in the function defineglobalvarnames:

    def defineinputdatavarnames(self):

        return {'SODA3': ['temp', 'salt', 'ssh', 'u', 'v','ui','vi','cn','hi','hs'],
                'GLORYS': ['thetao', 'so', 'zos', 'uo', 'vo', 'usi', 'vsi', 'siconc', 'sithick']}[self.indatatype]
        #return {'SODA3': ['temp', 'salt', 'ssh', 'u', 'v','ui','vi','cn','hi','hs'],            
        #        'GLORYS': ['temperature', 'salinity', 'ssh', 'u', 'v']}[self.indatatype]

    def defineromsgridpath(self):
        return {'antarctic': '/home/ueslei/Documents/model2roms/grid/antarctic.nc'}[self.outgrid]

    def defineforcingdatapath(self):
        return {'SODA3': "/home/ueslei/Documents/model2roms/input/",
                'GLORYS': "/home/ueslei/Documents/model2roms/input/"}[self.indatatype]

    def __init__(self):
        print('model2roms started ' + time.ctime(time.time()))
        os.environ['WRAP_STDERR'] = 'true'

        # Set compileAll to True if you want automatic re-compilation of all the
        # fortran files necessary to run model2roms. Options are "gfortran" or "ifort". Edit
        # compile.py to add other Fortran compilers.
        self.compileall = False
         # Create the bry, init, and clim files for a given grid and input data
        self.createoceanforcing = True
        # Create atmospheric forcing for the given grid
        self.createatmosforcing = False  # currently in beta stages and unavailable
        # Write ice values to file (for Arctic regions)
        self.writeice = True
        # ROMS sometimes requires input of ice and ssh, but if you dont have these write zero files to file
        self.set2DvarsToZero= False
        # Use ESMF for the interpolation. This requires that you have ESMF and ESMPy installed (import ESMF)
        self.useesmf = True
        # Apply filter to smooth the 2D fields after interpolation (time consuming but enhances results)
        self.usefilter = False
        # Format to write the ouput to: 'NETCDF4', 'NETCDF4_CLASSIC', 'NETCDF3_64BIT', or 'NETCDF3_CLASSIC'
        # Using NETCDF4 automatically turns on compression of files (ZLIB)
        self.myformat = 'NETCDF4'
        self.myzlib = True
        # Frequency of the input data: usually monthly
        self.timefrequencyofinputdata = "day"  # , "month", "hour"
        #  Define what grid type you wnat to interpolate from (input MODEL data)
        # Options:
        # 1. SODA3, GLORYS
        self.indatatype = 'GLORYS'
        # Define the names of the geographical variables in the input files
        if self.indatatype == "SODA3":
            self.grdtype     = 'regular'
            self.lonname     = "xt_ocean"
            self.latname     = "yt_ocean"
            self.depthname   = "st_ocean"
            self.lonname_u   = "xu_ocean"
            self.latname_u   = "yu_ocean"
            self.lonname_v   = "xu_ocean"
            self.latname_v   = "yu_ocean"
            self.timename    = "time"
            self.realm       = "ocean"
            self.fillvaluein = -1.e20
        if self.indatatype == "GLORYS":   
            self.grdtype     = 'regular'
            self.lonname     = "longitude"
            self.latname     = "latitude"
            self.depthname   = "depth"
            self.lonname_u   = "longitude"
            self.latname_u   = "latitude"
            self.lonname_v   = "longitude"
            self.latname_v   = "latitude"
            self.timename    = "time"
            self.realm       = "ocean"
            self.fillvaluein = -1.e20
        # Define contact info for final NetCDF files
        self.authorname  = "Ueslei Adriano Sutil"
        self.authoremail = "ueslei.sutil (at) inpe.br"
        # Define what grid type you want to interpolate from: Can be Z for SIGMA for ROMS
        # vertical coordinate system or ZLEVEL. Also define the name of the dimensions in the input files.
        # Options:
        # 1. SIGMA (not propoerly implemented yet), 2. ZLEVEL
        self.ingridtype = "ZLEVEL"
        # Define what grid type you want to interpolate to
        # Options: This is just the name of your grid used to identify your selection later
        self.outgrid     = "antarctic"
        self.outgridtype = "ROMS"
        # Define nmber of output depth levels
        self.nlevels     = 30
        # Define the grid stretching properties (leave default if uncertain what to pick)
        self.vstretching = 4
        self.vtransform  = 2
        self.theta_s     = 7.0
        self.theta_b     = 0.2
        self.tcline      = 20
        self.hc          = 20
        # Define the period to create forcing for
        self.start_year  = 2021
        self.end_year    = 2021
        self.start_month = 1
        self.end_month   = 1
        self.start_day   = 17
        self.end_day     = 18 

        if int(calendar.monthrange(self.start_year, self.start_month)[1]) < self.start_day:
            self.start_day = int(calendar.monthrange(self.start_year, self.start_month)[1])

        if int(calendar.monthrange(self.end_year, self.end_month)[1]) < self.end_day:
            self.end_day = int(calendar.monthrange(self.end_year, self.end_month)[1])

        self.startdate = datetime(self.start_year, self.start_month, self.start_day)
        self.enddate   = datetime(self.end_year, self.end_month, self.end_day)
        self.years     = [self.start_year + year for year in range(self.end_year + 1 - self.start_year)]
        
        self.globalvarnames    = self.defineglobalvarnames()
        self.inputdatavarnames = self.defineinputdatavarnames()
    
        self.modelpath    = self.defineforcingdatapath()
        self.romsgridpath = self.defineromsgridpath()

        if self.compileall is True:
            import compile;
            compile.compileallgfortran()

        if self.createatmosforcing or self.createoceanforcing:
            self.abbreviation = self.defineabbreviation()

            self.climname, self.initname, self.bryname = self.defineoutputfilenames()

            self.showinfo()

            if self.useesmf:
                try:
                    import ESMF
                except ImportError:
                    raise ImportError("Unable to import ESMF")
                ESMF.Manager(debug=True)

            # Create the grid object for the output grid
            self.grdROMS             = grd.Grd("ROMS", self)
            self.grdROMS.nlevels     = self.nlevels
            self.grdROMS.vstretching = self.vstretching
            self.grdROMS.vtransform  = self.vtransform
            self.grdROMS.theta_s     = self.theta_s
            self.grdROMS.theta_b     = self.theta_b
            self.grdROMS.tcline      = self.tcline
            self.grdROMS.hc          = self.hc
            self.grdROMS.lonname     = 'lon_rho'
            self.grdROMS.latname     = 'lat_rho'

            self.grdROMS.opennetcdf(self.romsgridpath)
            self.grdROMS.createobject(self)
            self.grdROMS.getdims()

            # Create the grid object for the input grid
            self.grdMODEL           = grd.Grd("FORCINGDATA", self)
            self.grdMODEL.grdType   = self.grdtype
            self.grdMODEL.lonName   = self.lonname
            self.grdMODEL.latName   = self.latname
            self.grdMODEL.depthName = self.depthname
            self.grdMODEL.fillval   = self.fillvaluein
