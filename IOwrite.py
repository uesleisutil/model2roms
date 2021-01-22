from datetime import datetime
from netCDF4 import Dataset
from netCDF4 import num2date
import numpy as np
import time
import os

"""
Created by Trond Kristiansen
https://github.com/trondkr/model2roms
"""

def help():
    """
    This function generates a CLIM file from scratch. The variables created include:
    salt, temp, u, v, ubar, vbar, zeta, and time. Time dimension for each variable is ocean_time which is days
    since 1948/1/1.

    """


def writeclimfile(confM2R, ntime, myvar, data1=None, data2=None, data3=None, data4=None):
    if confM2R.myformat == 'NETCDF4':
        myzlib = True
    else:
        myzlib = False

    grdROMS = confM2R.grdROMS

    if confM2R.grdROMS.ioClimInitialized is False:
        confM2R.grdROMS.ioClimInitialized = True
        if os.path.exists(confM2R.climname):
            os.remove(confM2R.climname)

        f1             = Dataset(confM2R.climname, mode='w', format=confM2R.myformat)
        f1.title       = "Climatology forcing file (CLI)"
        f1.description = "Created for grid file: %s" % (confM2R.abbreviation)
        f1.grd_file    = "Gridfile: %s" % (confM2R.romsgridpath)
        f1.history     = "Created " + time.ctime(time.time())
        f1.source      = "{} ({})".format(confM2R.authorname,confM2R.authoremail)

        # Define dimensions
        f1.createDimension('xi_rho', grdROMS.xi_rho)
        f1.createDimension('eta_rho', grdROMS.eta_rho)
        f1.createDimension('xi_u', grdROMS.xi_u)
        f1.createDimension('eta_u', grdROMS.eta_u)
        f1.createDimension('xi_v', grdROMS.xi_v)
        f1.createDimension('eta_v', grdROMS.eta_v)
        f1.createDimension('xi_psi', grdROMS.xi_psi)
        f1.createDimension('eta_psi', grdROMS.eta_psi)
        f1.createDimension('s_rho', len(grdROMS.s_rho))
        f1.createDimension('s_w', len(grdROMS.s_w))
        f1.createDimension('ocean_time', None)

        vnc               = f1.createVariable('lon_rho', 'd', ('eta_rho', 'xi_rho',), zlib=myzlib, fill_value=grdROMS.fillval)
        vnc.long_name     = 'Longitude of RHO-points'
        vnc.units         = 'degree_east'
        vnc.standard_name = 'longitude'
        vnc[:, :]         = grdROMS.lon_rho

        vnc               = f1.createVariable('lat_rho', 'd', ('eta_rho', 'xi_rho',), zlib=myzlib, fill_value=grdROMS.fillval)
        vnc.long_name     = 'Latitude of RHO-points'
        vnc.units         = 'degree_north'
        vnc.standard_name = 'latitude'
        vnc[:, :]         = grdROMS.lat_rho

        vnc               = f1.createVariable('lon_u', 'd', ('eta_u', 'xi_u',), zlib=myzlib, fill_value=grdROMS.fillval)
        vnc.long_name     = 'Longitude of U-points'
        vnc.units         = 'degree_east'
        vnc.standard_name = 'longitude'
        vnc[:, :]         = grdROMS.lon_u

        vnc               = f1.createVariable('lat_u', 'd', ('eta_u', 'xi_u',), zlib=myzlib, fill_value=grdROMS.fillval)
        vnc.long_name     = 'Latitude of U-points'
        vnc.units         = 'degree_north'
        vnc.standard_name = 'latitude'
        vnc[:, :]         = grdROMS.lat_u

        vnc               = f1.createVariable('lon_v', 'd', ('eta_v', 'xi_v',), zlib=myzlib, fill_value=grdROMS.fillval)
        vnc.long_name     = 'Longitude of V-points'
        vnc.units         = 'degree_east'
        vnc.standard_name = 'longitude'
        vnc[:, :]         = grdROMS.lon_v

        vnc               = f1.createVariable('lat_v', 'd', ('eta_v', 'xi_v',), zlib=myzlib, fill_value=grdROMS.fillval)
        vnc.long_name     = 'Latitude of V-points'
        vnc.units         = 'degree_north'
        vnc.standard_name = 'latitude'
        vnc[:, :]         = grdROMS.lat_v

        vnc               = f1.createVariable('lat_psi', 'd', ('eta_psi', 'xi_psi',), zlib=myzlib, fill_value=grdROMS.fillval)
        vnc.long_name     = 'Latitude of PSI-points'
        vnc.units         = 'degree_north'
        vnc.standard_name = 'latitude'
        vnc[:, :]         = grdROMS.lat_psi

        vnc               = f1.createVariable('lon_psi', 'd', ('eta_psi', 'xi_psi',), zlib=myzlib, fill_value=grdROMS.fillval)
        vnc.long_name     = 'Longitude of PSI-points'
        vnc.units         = 'degree_east'
        vnc.standard_name = 'longitude'
        vnc[:, :]         = grdROMS.lon_psi

        vnc           = f1.createVariable('h', 'd', ('eta_rho', 'xi_rho',), zlib=myzlib, fill_value=grdROMS.fillval)
        vnc.long_name = 'Bathymetry at RHO-points'
        vnc.units     = 'meter'
        vnc.field     = "bath, scalar"
        vnc[:, :]     = grdROMS.h

        vnc           = f1.createVariable('f', 'd', ('eta_rho', 'xi_rho',), zlib=myzlib, fill_value=grdROMS.fillval)
        vnc.long_name = 'Coriolis parameter at RHO-points'
        vnc.units     = 'second-1'
        vnc.field     = "Coriolis, scalar"
        vnc[:, :]     = grdROMS.f

        vnc           = f1.createVariable('pm', 'd', ('eta_rho', 'xi_rho',), zlib=myzlib, fill_value=grdROMS.fillval)
        vnc.long_name = 'curvilinear coordinate metric in XI'
        vnc.units     = 'meter-1'
        vnc.field     = "pm, scalar"
        vnc[:, :]     = grdROMS.pm

        vnc           = f1.createVariable('pn', 'd', ('eta_rho', 'xi_rho',), zlib=myzlib, fill_value=grdROMS.fillval)
        vnc.long_name = 'curvilinear coordinate metric in ETA'
        vnc.units     = 'meter-1'
        vnc.field     = "pn, scalar"
        vnc[:, :]     = grdROMS.pn

        vnc           = f1.createVariable('s_rho', 'd', ('s_rho',), zlib=myzlib, fill_value=grdROMS.fillval)
        vnc.long_name = "S-coordinate at RHO-points"
        vnc.valid_min = -1.
        vnc.valid_max = 0.
        if grdROMS.vtransform == 2:
            vnc.standard_name = "ocean_s_coordinate_g2"
            vnc.formula_terms = "s: s_rho C: Cs_r eta: zeta depth: h depth_c: hc"
        if grdROMS.vtransform == 1:
            vnc.standard_name = "ocean_s_coordinate_g1"
            vnc.formula_terms = "s: s_rho C: Cs_r eta: zeta depth: h depth_c: hc"
        vnc.field    = "s_rho, scalar"
        vnc[:]       = grdROMS.s_rho

        vnc           = f1.createVariable('s_w', 'd', ('s_w',), zlib=myzlib, fill_value=grdROMS.fillval)
        vnc.long_name = "S-coordinate at W-points"
        vnc.valid_min = -1.
        vnc.valid_max = 0.
        if grdROMS.vtransform == 2:
            vnc.standard_name = "ocean_s_coordinate_g2"
            vnc.formula_terms = "s: s_w C: Cs_w eta: zeta depth: h depth_c: hc"
        if grdROMS.vtransform == 1:
            vnc.standard_name = "ocean_s_coordinate_g1"
            vnc.formula_terms = "s: s_w C: Cs_w eta: zeta depth: h depth_c: hc"
        vnc.field     = "s_w, scalar"
        vnc[:]        = grdROMS.s_w

        vnc           = f1.createVariable('Cs_r', 'd', ('s_rho',), zlib=myzlib, fill_value=grdROMS.fillval)
        vnc.long_name = "S-coordinate stretching curves at RHO-points"
        vnc.valid_min = -1.
        vnc.valid_max = 0.
        vnc.field     = "s_rho, scalar"
        vnc[:]        = grdROMS.Cs_rho

        vnc           = f1.createVariable('Cs_w', 'd', ('s_w',), zlib=myzlib, fill_value=grdROMS.fillval)
        vnc.long_name = "S-coordinate stretching curves at W-points"
        vnc.valid_min = -1.
        vnc.valid_max = 0.
        vnc.field     = "s_w, scalar"
        vnc[:]        = grdROMS.Cs_w

        vnc           = f1.createVariable('hc', 'd')
        vnc.long_name = "S-coordinate parameter, critical depth";
        vnc.units     = "meter"
        vnc[:]        = grdROMS.hc

        vnc           = f1.createVariable('z_r', 'd', ('s_rho', 'eta_rho', 'xi_rho',), zlib=myzlib, fill_value=grdROMS.fillval)
        vnc.long_name = "Sigma layer to depth matrix";
        vnc.units     = "meter"
        vnc[:, :, :]  = grdROMS.z_r

        vnc           = f1.createVariable('z_w', 'd', ('s_w', 'eta_rho', 'xi_rho',), zlib=myzlib, fill_value=grdROMS.fillval)
        vnc.long_name = "Sigma layer to depth matrix";
        vnc.units     = "meter"
        vnc[:, :, :]  = grdROMS.z_w

        vnc           = f1.createVariable('Tcline', 'd')
        vnc.long_name = "S-coordinate surface/bottom layer width"
        vnc.units     = "meter"
        vnc[:]        = grdROMS.tcline

        vnc           = f1.createVariable('theta_s', 'd')
        vnc.long_name = "S-coordinate surface control parameter"
        vnc[:]        = grdROMS.theta_s

        vnc           = f1.createVariable('theta_b', 'd')
        vnc.long_name = "S-coordinate bottom control parameter"
        vnc[:]        = grdROMS.theta_b

        vnc           = f1.createVariable('angle', 'd', ('eta_rho', 'xi_rho',), zlib=myzlib, fill_value=grdROMS.fillval)
        vnc.long_name = "angle between xi axis and east"
        vnc.units     = "radian"
        vnc[:, :]     = grdROMS.angle

        # Now start creating variables for regular climatology/bry/init creations

        v_time           = f1.createVariable('ocean_time', 'd', ('ocean_time',), zlib=myzlib, fill_value=grdROMS.fillval)
        v_time.long_name = 'days since 1948-01-01 00:00:00'
        v_time.units     = 'days since 1948-01-01 00:00:00'
        v_time.field     = 'time, scalar, series'
        v_time.calendar  = 'standard'

        v_u           = f1.createVariable('u', 'f', ('ocean_time', 's_rho', 'eta_u', 'xi_u',), zlib=myzlib,fill_value=grdROMS.fillval)
        v_u.long_name = "u-momentum component"
        v_u.units     = "meter second-1"
        v_u.time      = "ocean_time"
        v_u.field     = "u-velocity, scalar, series"

        v_v           = f1.createVariable('v', 'f', ('ocean_time', 's_rho', 'eta_v', 'xi_v',), zlib=myzlib,fill_value=grdROMS.fillval)
        v_v.long_name = "v-momentum component"
        v_v.units     = "meter second-1"
        v_v.time      = "ocean_time"
        v_v.field     = "v-velocity, scalar, series"

        v_salt           = f1.createVariable('salt', 'f', ('ocean_time', 's_rho', 'eta_rho', 'xi_rho',), zlib=myzlib,fill_value=grdROMS.fillval)
        v_salt.long_name = "salinity"
        v_salt.time      = "ocean_time"
        v_salt.field     = "salinity, scalar, series"

        v_temp           = f1.createVariable('temp', 'f', ('ocean_time', 's_rho', 'eta_rho', 'xi_rho',), zlib=myzlib,fill_value=grdROMS.fillval)
        v_temp.long_name = "potential temperature"
        v_temp.units     = "Celsius"
        v_temp.time      = "ocean_time"
        v_temp.field     = "temperature, scalar, series"

        v_ssh           = f1.createVariable('zeta', 'f', ('ocean_time', 'eta_rho', 'xi_rho',), zlib=myzlib,fill_value=grdROMS.fillval)
        v_ssh.long_name = "sea level"
        v_ssh.units     = "meter"
        v_ssh.time      = "ocean_time"
        v_ssh.field     = "sea level, scalar, series"
  

        v_ubar           = f1.createVariable('ubar', 'f', ('ocean_time', 'eta_u', 'xi_u',), zlib=myzlib,fill_value=grdROMS.fillval)
        v_ubar.long_name = "u-2D momentum"
        v_ubar.units     = "meter second-1"
        v_ubar.time      = "ocean_time"
        v_ubar.field     = "u2-D velocity, scalar, series"

        v_vbar           = f1.createVariable('vbar', 'f', ('ocean_time', 'eta_v', 'xi_v',), zlib=myzlib,fill_value=grdROMS.fillval)
        v_vbar.long_name = "v-2D momentum"
        v_vbar.units     = "meter second-1"
        v_vbar.time      = "ocean_time"
        v_vbar.field     = "v2-D velocity, scalar, series"

        if confM2R.writeice:
            ageice           = f1.createVariable('ageice', 'f', ('ocean_time', 'eta_rho', 'xi_rho',), zlib=myzlib,fill_value=grdROMS.fillval)
            ageice.long_name = "time-averaged age of the ice"
            ageice.units     = "years"
            ageice.time      = "ocean_time"
            ageice.field     = "ice age, scalar, series"

            uice           = f1.createVariable('uice', 'd', ('ocean_time', 'eta_u', 'xi_u',), zlib=myzlib,fill_value=grdROMS.fillval)
            uice.long_name = "time-averaged u-component of ice velocity"
            uice.units     = "meter second-1"
            uice.time      = "ocean_time"
            uice.field     = "u-component of ice velocity, scalar, series"

            vice           = f1.createVariable('vice', 'd', ('ocean_time', 'eta_v', 'xi_v',), zlib=myzlib,fill_value=grdROMS.fillval)
            vice.long_name = "time-averaged v-component of ice velocity"
            vice.units     = "meter second-1"
            vice.time      = "ocean_time"
            vice.field     = "v-component of ice velocity, scalar, series"

            aice           = f1.createVariable('aice', 'f', ('ocean_time', 'eta_rho', 'xi_rho',), zlib=myzlib,fill_value=grdROMS.fillval)
            aice.long_name = "time-averaged fraction of cell covered by ice"
            aice.time      = "ocean_time"
            aice.field     = "ice concentration, scalar, series"

            hice           = f1.createVariable('hice', 'f', ('ocean_time', 'eta_rho', 'xi_rho',), zlib=myzlib,fill_value=grdROMS.fillval)
            hice.long_name = "time-averaged average ice thickness in cell"
            hice.units     = "meter"
            hice.time      = "ocean_time"
            hice.field     = "ice thickness, scalar, series"

            snow_thick           = f1.createVariable('snow_thick', 'f', ('ocean_time', 'eta_rho', 'xi_rho',), zlib=myzlib,fill_value=grdROMS.fillval)
            snow_thick.long_name = "time-averaged thickness of snow cover"
            snow_thick.units     = "meter"
            snow_thick.time      = "ocean_time"
            snow_thick.field     = "snow thickness, scalar, series"

            ti           = f1.createVariable('ti', 'f', ('ocean_time', 'eta_rho', 'xi_rho',), zlib=myzlib,fill_value=grdROMS.fillval)
            ti.long_name = "time-averaged interior ice temperature"
            ti.units     = "degrees Celcius"
            ti.time      = "ocean_time"
            ti.field     = "interior temperature, scalar, series"

            sfwat           = f1.createVariable('sfwat', 'f', ('ocean_time', 'eta_rho', 'xi_rho',), zlib=myzlib,fill_value=grdROMS.fillval)
            sfwat.long_name = "time-averaged surface melt water thickness on ice"
            sfwat.units     = "meter"
            sfwat.time      = "ocean_time"
            sfwat.field     = "melt water thickness, scalar, series"

            tisrf           = f1.createVariable('tisrf', 'f', ('ocean_time', 'eta_rho', 'xi_rho',), zlib=myzlib,fill_value=grdROMS.fillval)
            tisrf.long_name = "time-averaged temperature of ice surface"
            tisrf.units     = "degrees Celcius"
            tisrf.time      = "ocean_time"
            tisrf.field     = "surface temperature, scalar, series"

            sig11           = f1.createVariable('sig11', 'f', ('ocean_time', 'eta_rho', 'xi_rho',), zlib=myzlib,fill_value=grdROMS.fillval)
            sig11.long_name = "time-averaged internal ice stress 11 component"
            sig11.units     = "Newton meter-1"
            sig11.time      = "ocean_time"
            sig11.field     = "ice stress 11, scalar, series"

            sig12           = f1.createVariable('sig12', 'f', ('ocean_time', 'eta_rho', 'xi_rho',), zlib=myzlib,fill_value=grdROMS.fillval)
            sig12.long_name = "time-averaged internal ice stress 12 component"
            sig12.units     = "Newton meter-1"
            sig12.time      = "ocean_time"
            sig12.field     = "ice stress 12, scalar, series"

            sig22           = f1.createVariable('sig22', 'f', ('ocean_time', 'eta_rho', 'xi_rho',), zlib=myzlib,fill_value=grdROMS.fillval)
            sig22.long_name = "time-averaged internal ice stress 22 component"
            sig22.units     = "Newton meter-1"
            sig22.time      = "ocean_time"
            sig22.field     = "ice stress 22, scalar, series"

    else:
        f1 = Dataset(confM2R.climname, mode='a', format=confM2R.myformat)

    if myvar == confM2R.globalvarnames[0]:

        if grdROMS.timeunits[0:7] == "seconds":
            f1.variables['ocean_time'][ntime] = grdROMS.time
            d = num2date(grdROMS.time, units=f1.variables['ocean_time'].long_name,
                            calendar=f1.variables['ocean_time'].calendar)
        else:
            f1.variables['ocean_time'][ntime] = grdROMS.time

            d = num2date(grdROMS.time, units=f1.variables['ocean_time'].long_name,
                            calendar=f1.variables['ocean_time'].calendar)
        grdROMS.message = d

    if myvar == 'temperature':
        f1.variables['temp'][ntime, :, :, :] = data1
    if myvar == 'salinity':
        f1.variables['salt'][ntime, :, :, :] = data1
    if myvar == 'ssh':
        f1.variables['zeta'][ntime, :, :] = data1
    if myvar == 'vvel':
        f1.variables['u'][ntime, :, :, :] = data1
        f1.variables['v'][ntime, :, :, :] = data2

        f1.variables['ubar'][ntime, :, :] = data3
        f1.variables['vbar'][ntime, :, :] = data4

    if confM2R.writeice:
        if myvar == "ageice":
            data1 = np.where(abs(data1) > 100, 0, data1)
            f1.variables['ageice'][ntime, :, :] = data1
            #f1.variables['ageice'][ntime, :, :]     = 0.

        if myvar == 'uice':
            data1 = np.where(abs(data1) > 120, 0, data1)
            f1.variables['uice'][ntime, :, :]  = data1
            f1.variables['sfwat'][ntime, :, :] = 0.
            f1.variables['tisrf'][ntime, :, :] = 0.
            f1.variables['ti'][ntime, :, :]    = 0.
            f1.variables['sig11'][ntime, :, :] = 0.
            f1.variables['sig12'][ntime, :, :] = 0.
            f1.variables['sig22'][ntime, :, :] = 0.

           # if confM2R.indatatype == 'GLORYS':
                # Special care for GLORYS as dataset does not contain sea ice age and snow thickness
                #f1.variables['ageice'][ntime, :, :]     = 0.
                #f1.variables['snow_thick'][ntime, :, :] = 0

        if myvar == 'vice':
            data1 = np.where(abs(data1) > 120, 0, data1)
            f1.variables['vice'][ntime, :, :] = data1
        if myvar == 'aice':
            data1 = np.where(abs(data1) > 120, 0, data1)
            f1.variables['aice'][ntime, :, :] = data1
        if myvar == 'hice':
            data1 = np.where(abs(data1) > 10, 0, data1)
            #data1 = np.ma.masked_where(abs(data1) > 10, data1)
            f1.variables['hice'][ntime, :, :] = data1
        if myvar == 'snow_thick':
            #data1 = np.ma.masked_where(abs(data1) > 100, data1)
            data1 = np.where(abs(data1) > 10, 0, data1)
            f1.variables['snow_thick'][ntime, :, :] = data1

    f1.close()
