from datetime import datetime
from netCDF4 import Dataset
import numpy as np

import IOverticalGrid

try:
    import ESMF
except ImportError:
    print("Could not find module ESMF")
    pass

"""
Trond Kristiansen
https://github.com/trondkr/model2roms
"""

class Grd:

    def __init__(self, grdtype, confM2R):
        """
        The object is initialised and created through the __init__ method
        As an example of how to use, these lines return a grid object called grdTEST:
        => import grd
        => grdTEST = grd.grdClass("grdfilename","ROMS")
        """
        self.type = grdtype
        self.grdName = confM2R.outgrid
        self.realm = confM2R.realm
        self.grdfilename = None

    def opennetcdf(self, grdfilename):
        
        self.grdfilename = grdfilename
        """Open the netCDF file and store the contents in arrays associated with variable names"""
        try:
            self.cdf = Dataset(self.grdfilename, "r")

        except IOError:
            print('Could not open file {}'.format(self.grdfilename))
            print('Exception caught in: opennetcdf(grdfilename)')

    def createobject(self, confM2R):
        """
        This method creates a new object by reading the grd input file. All
        dimensions (eta, xi, lon, lat etc.) are defined here and used througout these scripts.
        Also, the depth matrix is calculated in this function by calling IOverticalGrid.py (ROMS grid only). For
        input model depths, the depth array is a one dimensional. If input data has 2 or 3 dimensions, this
        has to be accounted for througout the soda2roms package as one dimension is currently only supported.
        """
        if self.type == 'FORCINGDATA':
            self.lon = self.cdf.variables[str(confM2R.lonname)][:]
            self.lat = self.cdf.variables[str(confM2R.latname)][:]
            self.h = self.cdf.variables[str(confM2R.depthname)][:]
            self.nlevels = len(self.h)
            self.fillval = -9.99e+33

            if self.lon.ndim == 1:
                self.lon, self.lat = np.meshgrid(self.lon, self.lat)

            # Create grid for ESMF interpolation
            if confM2R.useesmf:
                self.esmfgrid = ESMF.Grid(filename=self.grdfilename, filetype=ESMF.FileFormat.GRIDSPEC,
                                          is_sphere=True, coord_names=[str(confM2R.lonname), str(confM2R.latname)],
                                          add_mask=False)
                self.esmfgrid_u = ESMF.Grid(filename=self.grdfilename, filetype=ESMF.FileFormat.GRIDSPEC,
                                          is_sphere=True, coord_names=[str(confM2R.lonname_u), str(confM2R.latname_u)],
                                          add_mask=False)
                self.esmfgrid_v = ESMF.Grid(filename=self.grdfilename, filetype=ESMF.FileFormat.GRIDSPEC,
                                          is_sphere=True, coord_names=[str(confM2R.lonname_v), str(confM2R.latname_v)],
                                          add_mask=False)

            if confM2R.indatatype == 'SODA3':
                self.fillval = -1.e+20
            if confM2R.indatatype == 'GLORYS':
                self.fillval = 9.96921e+36
                   
            IOverticalGrid.get_z_levels(self)

        if self.type in ['ROMS']:

            self.write_clim     = True
            self.write_bry      = True
            self.write_init     = True

            self.lonname = 'lon_rho'
            self.latname = 'lat_rho'

            """Set initTime to 1 if you dont want the first timestep to be
            the initial field (no ubar and vbar if time=0)"""

            self.inittime    = 0
            self.ocean_time  = 0
            self.NT          = 2
            self.tracer      = self.NT

            self.message     = None
            self.time        = 0
            self.reftime     = 0
            self.grdtype     = 'regular'
            self.mask_rho    = self.cdf.variables["mask_rho"][:, :]
            self.lon_rho     = self.cdf.variables["lon_rho"][:, :]
            self.lat_rho     = self.cdf.variables["lat_rho"][:, :]
            self.h           = self.cdf.variables["h"][:, :]
            self.hmin        = self.h[self.h > 0].min()
            self.vtransform  = confM2R.vtransform
            self.nlevels     = confM2R.nlevels
            self.vstretching = confM2R.vstretching     
            self.theta_s     = confM2R.theta_s
            self.theta_b     = confM2R.theta_b
            self.tcline      = confM2R.tcline
            self.hc          = confM2R.hc
            if self.vtransform == 1:
                self.hc = min(self.hmin, self.tcline)
                self.hc = self.tcline
                if (self.tcline > self.hmin):
                    print('Vertical transformation parameters are not defined correctly in either gridid.txt or in the history files: \n Tc\
line = %d and hmin = %d. \n You need to make sure that tcline <= hmin when using transformation 1.' % (
                    self.tcline, self.hmin))
            else:
                self.hc = self.tcline

            zeta = None
            if zeta is None:
                self.zeta = np.zeros(self.h.shape)
            else:
                self.zeta = zeta

            # for findvar in self.cdf.variables:
            #    if findvar=="hraw":
            #        self.hraw     = self.cdf.variables["hraw"][:,:,:]

            self.lon_u = self.cdf.variables["lon_u"][:, :]
            self.lat_u = self.cdf.variables["lat_u"][:, :]
            self.mask_u = self.cdf.variables["mask_u"][:, :]
            for findvar in self.cdf.variables:
                if findvar == "lon_vert":
                    self.lon_vert = self.cdf.variables["lon_vert"][:, :]
                    self.lat_vert = self.cdf.variables["lat_vert"][:, :]

            for findvar in self.cdf.variables:
                if findvar == "x_rho":
                    self.x_rho = self.cdf.variables["x_rho"][:, :]
                    self.y_rho = self.cdf.variables["y_rho"][:, :]

            for findvar in self.cdf.variables:
                if findvar == "x_u":
                    self.x_u = self.cdf.variables["x_u"][:, :]
                    self.y_u = self.cdf.variables["y_u"][:, :]

            for findvar in self.cdf.variables:
                if findvar == "x_v":
                    self.x_v = self.cdf.variables["x_v"][:, :]
                    self.y_v = self.cdf.variables["y_v"][:, :]

            for findvar in self.cdf.variables:
                if findvar == "x_psi":
                    self.x_psi = self.cdf.variables["x_psi"][:, :]
                    self.y_psi = self.cdf.variables["y_psi"][:, :]

            for findvar in self.cdf.variables:
                if findvar == "x_vert":
                    self.x_vert = self.cdf.variables["x_vert"][:, :]
                    self.y_vert = self.cdf.variables["y_vert"][:, :]

            for findvar in self.cdf.variables:
                if findvar == "xl":
                    self.xl = self.cdf.variables["xl"][:]
                    self.el = self.cdf.variables["el"][:]

            for findvar in self.cdf.variables:
                if findvar == "dmde":
                    self.dmde = self.cdf.variables["dmde"][:, :]
                    self.dndx = self.cdf.variables["dndx"][:, :]

            self.lon_v = self.cdf.variables["lon_v"][:, :]
            self.lat_v = self.cdf.variables["lat_v"][:, :]
            self.mask_v = self.cdf.variables["mask_v"][:, :]

            self.spherical = self.cdf.variables["spherical"][:]

            self.lon_psi = self.lon_u[:-1, :]
            self.lat_psi = self.lat_v[:, :-1]
            self.mask_psi = self.mask_v[:, :-1]

            self.f = self.cdf.variables["f"][:, :]
            self.angle = self.cdf.variables["angle"][:, :]

            self.pm = self.cdf.variables["pm"][:, :]
            self.invpm = 1.0 / np.asarray(self.cdf.variables["pm"][:, :])
            self.pn = self.cdf.variables["pn"][:, :]
            self.invpn = 1.0 / np.asarray(self.cdf.variables["pn"][:, :])

            self.Lp = len(self.lat_rho[1, :])
            self.Mp = len(self.lat_rho[:, 1])

            self.fillval = -9.99e33

            self.eta_rho = self.Mp
            self.eta_u = self.Mp
            self.eta_v = self.Mp - 1
            self.eta_psi = self.Mp - 1
            self.xi_rho = self.Lp
            self.xi_u = self.Lp - 1
            self.xi_v = self.Lp
            self.xi_psi = self.Lp - 1

            """Boolean to check if we need to initialize the CLIM file before writing"""
            self.ioClimInitialized = False
            self.ioInitInitialized = False

            if self.lon_rho.ndim == 1:
                self.lon_rho, self.lat_rho = np.meshgrid(self.lon_rho, self.lat_rho)
                self.lon_u, self.lat_u = np.meshgrid(self.lon_u, self.lat_u)
                self.lon_v, self.lat_v = np.meshgrid(self.lon_v, self.lat_v)

            """Setup the vertical coordinate system"""
            IOverticalGrid.calculateVgrid(self)

            if (confM2R.useesmf):
                self.esmfgrid_u = ESMF.Grid(filename=self.grdfilename, filetype=ESMF.FileFormat.GRIDSPEC,
                                            coord_names=['lon_u', 'lat_u'], add_mask=False)
                self.esmfgrid_v = ESMF.Grid(filename=self.grdfilename, filetype=ESMF.FileFormat.GRIDSPEC,
                                            is_sphere=True, coord_names=['lon_v', 'lat_v'], add_mask=False)
                self.esmfgrid = ESMF.Grid(filename=self.grdfilename, filetype=ESMF.FileFormat.GRIDSPEC,
                                          is_sphere=True, coord_names=[self.lonname, self.latname], add_mask=False)

    def getdims(self):
        if self.type in ["ROMS"]:
            self.Lp = len(self.lat_rho[1, :])
            self.Mp = len(self.lat_rho[:, 1])

        if self.type in ['FORCINGDATA']:
            self.Lp = len(self.lat[1, :])
            self.Mp = len(self.lat[:, 1])

        self.M = self.Mp - 1
        self.L = self.Lp - 1
