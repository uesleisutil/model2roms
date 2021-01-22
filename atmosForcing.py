import time
from datetime import datetime, timedelta
import os, sys, string
from netCDF4 import Dataset, num2date
import numpy as np
import IOatmos
import grd
import extrapolate as ex
try:
    import ESMF
except ImportError:
    print("Could not find module ESMF. Required")
    sys.exit()
"""
This funcion was created by Trond Kristiansen
https://github.com/trondkr/model2roms
"""
def laplaceFilter(field, threshold, toxi, toeta):
    undef = 2.0e+35
    tx    = 0.9*undef
    critx = 0.01
    cor   = 1.6
    mxs   = 10

    field = np.where(abs(field)>threshold,undef,field)

    field = ex.extrapolate.fill(int(1),int(toxi),
                                  int(1),int(toeta),
                                  float(tx), float(critx), float(cor), float(mxs),
                                  np.asarray(field, order='F'),
                                  int(toxi),
                                  int(toeta))
    return field

def createAtmosFileUV(grdROMS,modelpath,atmospath,startdate,enddate,useESMF,myformat,abbreviation,mytype,gridtype):

    # Setup 
    years = [(int(startdate.year) + kk) for kk in range(1 + int(enddate.year) - int(startdate.year))]
   
    # Create the objects for source and destination grids
   
    # Get the "Fraction of sfc area covered by ocean
    nor      = atmospath + "NRCP45AERCN_f19_g16_CLE_02.cam2.h0.2006-01.nc"
    cdf      = Dataset(nor,"r")
    OCENFRAC = cdf.variables["OCNFRAC"][:]
    cdf.close()
    Fill     = -999.0

    grdMODEL = grd.grdClass(nor, mytype, mytype, useESMF,'atmos')
    
    # Create the outputfile
    outfilename=  abbreviation + '_windUV_' + str(mytype) + '_' + str(startdate.year) + '_to_' + str(enddate.year) + '.nc'
    IOatmos.createNetCDFFileUV(grdROMS, outfilename, myformat, mytype)
    
    # Setup ESMF for interpolation (calculates weights)
    grdMODEL.fieldSrc = ESMF.Field(grdMODEL.esmfgrid, "fieldSrc", staggerloc=ESMF.StaggerLoc.CENTER)
    grdMODEL.fieldDst_rho = ESMF.Field(grdROMS.esmfgrid, "fieldDst", staggerloc=ESMF.StaggerLoc.CENTER)
    grdMODEL.regridSrc2Dst_rho = ESMF.Regrid(grdMODEL.fieldSrc, grdMODEL.fieldDst_rho, regrid_method=ESMF.RegridMethod.BILINEAR)
