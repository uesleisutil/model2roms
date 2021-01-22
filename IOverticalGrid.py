
from datetime import datetime
import numpy as np
import warnings

'''
Various vertical coordinates

Presently, only ocean s-coordinates are supported. Future plans will be to
include all of the vertical coordinate systems defined by the CF conventions.

vgrid.py function copied from https://github.com/kshedstrom/pyroms (Frederic Castruccio)

Adapted by Trond Kristiansen
https://github.com/trondkr/model2roms
'''

def calculateVgrid(self):
    if self.vtransform == 1:
        vgrid = s_coordinate(self.h, self.theta_b, self.theta_s, self.tcline, self.nlevels, self.vtransform, self.vstretching, zeta=None)
    elif self.vtransform == 2 and self.vstretching == 2:
        vgrid = s_coordinate_2(self.h, self.theta_b, self.theta_s, self.tcline, self.nlevels, self.vtransform, self.vstretching, zeta=None)
    elif self.vtransform == 2 and self.vstretching == 4:
        vgrid = s_coordinate_4(self.h, self.theta_b, self.theta_s, self.tcline, self.nlevels, self.vtransform, self.vstretching, zeta=None)
    else:
        raise Warning('Unknow vertical transformation Vtrans')
    self.z_r    = vgrid.z_r[0,:]
    self.z_w    = vgrid.z_w[0,:]
    self.Cs_rho = vgrid.Cs_r
    self.Cs_w   = vgrid.Cs_w
    self.s_rho  = vgrid.s_rho
    self.s_w    = vgrid.s_w

class s_coordinate(object):
    """
    Song and Haidvogel (1994) vertical coordinate transformation (Vtransform=1) and
    stretching functions (Vstretching=1).
    
    return an object that can be indexed to return depths

    s = s_coordinate(h, theta_b, theta_s, Tcline, N)
    """

    def __init__(self, h, theta_b, theta_s, tcline, N, vtransform, vstretching, zeta=None):
      
        self.h           = np.asarray(h)
        self.hmin        = h.min()
        self.theta_b     = theta_b
        self.theta_s     = theta_s
        self.tcline      = tcline
        self.N           = int(N)
        self.Np          = self.N+1
        self.vtransform  = vtransform
        self.vstretching = vstretching
        self.hc          = min(self.hmin, self.tcline)

        self.Vtrans = 1

        if self.vtransform==1:
            if (self.tcline > self.hmin):
                warnings.warn('Vertical transformation parameters are not defined correctly in either gridid.txt or in the history files: \n Tcline = %d and hmin = %d. \n You need to make sure that Tcline <= hmin when using transformation 1.' %(self.Tcline,self.hmin))

        self.c1 = 1.0
        self.c2 = 2.0
        self.p5 = 0.5       

        if zeta is None:
            self.zeta = np.zeros(h.shape)
        else:
            self.zeta = zeta
        
        self._get_s_rho()
        self._get_s_w()
        self._get_Cs_r()
        self._get_Cs_w()

        self.z_r = z_r(self.h, self.hc, self.N, self.s_rho, self.Cs_r, self.zeta, self.Vtrans)
        self.z_w = z_w(self.h, self.hc, self.Np, self.s_w, self.Cs_w, self.zeta, self.Vtrans)


    def _get_s_rho(self):
        lev        = np.arange(1,self.N+1,1)
        ds         = 1.0 / self.N
        self.s_rho = -self.c1 + (lev - self.p5) * ds

    def _get_s_w(self):
        lev = np.arange(0,self.Np,1)
        ds = 1.0 / (self.Np-1)
        self.s_w = -self.c1 + lev * ds

    def _get_Cs_r(self):
        if (self.theta_s >= 0):
            Ptheta = np.sinh(self.theta_s * self.s_rho) / np.sinh(self.theta_s)
            Rtheta = np.tanh(self.theta_s * (self.s_rho + self.p5)) / \
                      (self.c2 * np.tanh(self.p5 * self.theta_s)) - self.p5
            self.Cs_r = (self.c1 - self.theta_b) * Ptheta + self.theta_b * Rtheta
        else:
            self.Cs_r = self.s_rho

    def _get_Cs_w(self):
        if (self.theta_s >= 0):
            Ptheta = np.sinh(self.theta_s * self.s_w) / np.sinh(self.theta_s)
            Rtheta = np.tanh(self.theta_s * (self.s_w + self.p5)) / \
                      (self.c2 * np.tanh(self.p5 * self.theta_s)) - self.p5
            self.Cs_w = (self.c1 - self.theta_b) * Ptheta + self.theta_b * Rtheta
        else:
            self.Cs_w = self.s_w



class s_coordinate_2(s_coordinate):
    """
    A. Shchepetkin (2005) UCLA-ROMS vertical coordinate transformation (Vtransform=2) and
    stretching functions (Vstretching=2).
    
    return an object that can be indexed to return depths

    s = s_coordinate_2(h, theta_b, theta_s, Tcline, N)
    """

    def __init__(self, h, theta_b, theta_s, tcline, N,  vtransform, vstretching, zeta=None):
        self.h           = np.asarray(h)
        self.hmin        = h.min()
        self.theta_b     = theta_b
        self.theta_s     = theta_s
        self.tcline      = tcline
        self.N           = int(N)
        self.Np          = self.N+1
        self.vtransform  = vtransform
        self.vstretching = vstretching

        self.hc = self.tcline

        self.Vtrans = 2

        self.Aweight = 1.0
        self.Bweight = 1.0

        self.c1 = 1.0
        self.c2 = 2.0
        self.p5 = 0.5

        if zeta is None:
            self.zeta = np.zeros(h.shape)
        else:
            self.zeta = zeta

        self._get_s_rho()
        self._get_s_w()
        self._get_Cs_r()
        self._get_Cs_w()

        self.z_r = z_r(self.h, self.hc, self.N, self.s_rho, self.Cs_r, self.zeta, self.Vtrans)
        self.z_w = z_w(self.h, self.hc, self.Np, self.s_w, self.Cs_w, self.zeta, self.Vtrans)

    def _get_s_rho(self):
        super(s_coordinate_2, self)._get_s_rho()

    def _get_s_w(self):
        super(s_coordinate_2, self)._get_s_w()

    def _get_Cs_r(self):
        if (self.theta_s >= 0):
            Csur = (self.c1 - np.cosh(self.theta_s * self.s_rho)) / \
                     (np.cosh(self.theta_s) - self.c1)
            if (self.theta_b >= 0):
                Cbot = np.sinh(self.theta_b * (self.s_rho + self.c1)) / \
                       np.sinh(self.theta_b) - self.c1
                Cweight = (self.s_rho + self.c1)**self.Aweight * \
                          (self.c1 + (self.Aweight / self.Bweight) * \
                          (self.c1 - (self.s_rho + self.c1)**self.Bweight))
                self.Cs_r = Cweight * Csur + (self.c1 - Cweight) * Cbot
            else:
                self.Cs_r = Csur
        else:
            self.Cs_r = self.s_rho 
      
    def _get_Cs_w(self):
        if (self.theta_s >= 0):
            Csur = (self.c1 - np.cosh(self.theta_s * self.s_w)) / \
                     (np.cosh(self.theta_s) - self.c1)
            if (self.theta_b >= 0):
                Cbot = np.sinh(self.theta_b * (self.s_w + self.c1)) / \
                       np.sinh(self.theta_b) - self.c1
                Cweight = (self.s_w + self.c1)**self.Aweight * \
                          (self.c1 + (self.Aweight / self.Bweight) * \
                          (self.c1 - (self.s_w + self.c1)**self.Bweight))
                self.Cs_w = Cweight * Csur + (self.c1 - Cweight) * Cbot
            else:
                self.Cs_w = Csur
        else:
            self.Cs_w = self.s_w

class s_coordinate_4(s_coordinate):
    """
    A. Shchepetkin (2005) UCLA-ROMS vertical coordinate transformation (Vtransform=2) and
    stretching functions (Vstretching=4).
    
    return an object that can be indexed to return depths

    s = s_coordinate_4(h, theta_b, theta_s, Tcline, N)
    """

    def __init__(self, h, theta_b, theta_s, tcline, N,  vtransform, vstretching, zeta=None):
        self.h           = np.asarray(h)
        self.hmin        = h.min()
        self.theta_b     = theta_b
        self.theta_s     = theta_s
        self.tcline      = tcline
        self.N           = int(N)
        self.Np          = self.N+1
        self.vtransform  = vtransform
        self.vstretching = vstretching

        self.hc = self.tcline

        self.Vtrans = 4

        self.c1 = 1.0
        self.c2 = 2.0
        self.p5 = 0.5

        if zeta is None:
            self.zeta = np.zeros(h.shape)
        else:
            self.zeta = zeta

        self._get_s_rho()
        self._get_s_w()
        self._get_Cs_r()
        self._get_Cs_w()

        self.z_r = z_r(self.h, self.hc, self.N, self.s_rho, self.Cs_r, self.zeta, self.Vtrans)
        self.z_w = z_w(self.h, self.hc, self.Np, self.s_w, self.Cs_w, self.zeta, self.Vtrans)
        
    def _get_s_rho(self):
        super(s_coordinate_4, self)._get_s_rho()

    def _get_s_w(self):
        super(s_coordinate_4, self)._get_s_w()

    def _get_Cs_r(self):
        if (self.theta_s > 0):
            Csur = (self.c1 - np.cosh(self.theta_s * self.s_rho)) / \
                     (np.cosh(self.theta_s) - self.c1)
        else:
            Csur = -self.s_rho**2
        if (self.theta_b > 0):
            Cbot = (np.exp(self.theta_b * Csur) - self.c1 ) / \
                   (self.c1 - np.exp(-self.theta_b))
            self.Cs_r = Cbot
        else:
            self.Cs_r = Csur         

    def _get_Cs_w(self):
        if (self.theta_s > 0):
            Csur = (self.c1 - np.cosh(self.theta_s * self.s_w)) / \
                     (np.cosh(self.theta_s) - self.c1)
        else:
            Csur = -self.s_w**2
        if (self.theta_b > 0):
            Cbot = (np.exp(self.theta_b * Csur) - self.c1 ) / \
                   ( self.c1 - np.exp(-self.theta_b) )
            self.Cs_w = Cbot
        else:
            self.Cs_w = Csur


class z_r(object):
    """
    return an object that can be indexed to return depths of rho point

    z_r = z_r(h, hc, N, s_rho, Cs_r, zeta, Vtrans)
    """

    def __init__(self, h, hc, N, s_rho, Cs_r, zeta, Vtrans):
        self.h      = h
        self.hc     = hc
        self.N      = N
        self.s_rho  = s_rho
        self.Cs_r   = Cs_r
        self.zeta   = zeta
        self.Vtrans = Vtrans

    def __getitem__(self, key):

        if isinstance(key, tuple) and len(self.zeta.shape) > len(self.h.shape):
            zeta = self.zeta[key[0]]
            res_index = (slice(None),) + key[1:]
        elif len(self.zeta.shape) > len(self.h.shape):
            zeta = self.zeta[key]
            res_index = slice(None)
        else:
            zeta = self.zeta
            res_index = key
     
        if self.h.ndim == zeta.ndim:       # Assure a time-dimension exists
            zeta = zeta[np.newaxis, :]
        
        ti = zeta.shape[0]
        z_r = np.empty((ti, self.N) + self.h.shape, 'd')
        if self.Vtrans == 1:
            for n in range(ti):
                for  k in range(self.N):
                    z0 = self.hc * self.s_rho[k] + (self.h - self.hc) * self.Cs_r[k]
                    z_r[n,k,:] = z0 + zeta[n,:] * (1.0 + z0 / self.h)
        elif self.Vtrans == 2 or self.Vtrans == 4:
            for n in range(ti):
                for  k in range(self.N):
                    z0 = (self.hc * self.s_rho[k] + self.h * self.Cs_r[k]) / \
                          (self.hc + self.h)
                    z_r[n,k,:] = zeta[n,:] + (zeta[n,:] + self.h) * z0

        return np.squeeze(z_r[res_index])


class z_w(object):
    """
    return an object that can be indexed to return depths of w point

    z_w = z_w(h, hc, Np, s_w, Cs_w, zeta, Vtrans)
    """

    def __init__(self, h, hc, Np, s_w, Cs_w, zeta, Vtrans):
        self.h      = h
        self.hc     = hc
        self.Np     = Np
        self.s_w    = s_w
        self.Cs_w   = Cs_w
        self.zeta   = zeta
        self.Vtrans = Vtrans

    def __getitem__(self, key):

        if isinstance(key, tuple) and len(self.zeta.shape) > len(self.h.shape):
            zeta = self.zeta[key[0]]
            res_index = (slice(None),) + key[1:]
        elif len(self.zeta.shape) > len(self.h.shape):
            zeta = self.zeta[key]
            res_index = slice(None)
        else:
            zeta = self.zeta
            res_index = key
     
        if self.h.ndim == zeta.ndim:       # Assure a time-dimension exists
            zeta = zeta[np.newaxis, :]
      
        ti = zeta.shape[0]
        z_w = np.empty((ti, self.Np) + self.h.shape, 'd')
        if self.Vtrans == 1:
            for n in range(ti):
                for  k in range(self.Np):
                    z0 = self.hc * self.s_w[k] + (self.h - self.hc) * self.Cs_w[k]
                    z_w[n,k,:] = z0 + zeta[n,:] * (1.0 + z0 / self.h)
        elif self.Vtrans == 2 or self.Vtrans == 4:
            for n in range(ti):
                for  k in range(self.Np):
                    z0 = (self.hc * self.s_w[k] + self.h * self.Cs_w[k]) / \
                          (self.hc + self.h)
                    z_w[n,k,:] = zeta[n,:] + (zeta[n,:] + self.h) * z0

        return np.squeeze(z_w[res_index])



def get_z_levels(self):
    """
    Get a list of all the variables contained in netCDF file "fi32lename"
    """
    self.z_r=-self.h

    if len(self.z_r)==0:
        print(("No depth matrix found in file %s"%(self.selffilename)))


