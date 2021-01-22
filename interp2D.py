from __future__ import print_function
import numpy as np
import datetime
import extrapolate as ex
import os
try:
    import ESMF
except ImportError:
    print("Could not find module ESMF")
    pass

"""
Created by Trond Kristiansen
https://github.com/trondkr/model2roms
"""

def laplacefilter(field, threshold, toxi, toeta):
    undef = 2.0e+35
    tx    = 0.9 * undef
    critx = 0.01
    cor   = 1.6
    mxs   = 10

    field = np.where(abs(field) > threshold, undef, field)

    field = ex.extrapolate.fill(int(1), int(toxi),
                                int(1), int(toeta),
                                float(tx), float(critx), float(cor), float(mxs),
                                np.asarray(field, order='K'),
                                int(toxi),
                                int(toeta))
    return field


def dohorinterpolationregulargrid(confM2R, mydata, myvar):
    indexROMS_Z_ST, toxi, toeta, mymask = setupIndexes(confM2R, myvar)
    array1 = np.zeros((indexROMS_Z_ST), dtype=np.float)

    # 2D or 3D interpolation
    depthlevels = confM2R.grdMODEL.nlevels
    if myvar in ['ssh', 'ageice', 'uice', 'vice', 'aice', 'hice', 'snow_thick','hs']:
        depthlevels=1
    
    for k in range(depthlevels):
        if confM2R.useesmf:
            if depthlevels == 1:
                indata = np.squeeze(mydata[:, :])
            else:
                indata = np.squeeze(mydata[k, :, :])

            # We interpolate to RHO fields for all variables and then we later interpolate RHO points to U and V points
            # But input data are read on U and V and RHO grids if they differ..
            # EXCEPTION IS UICE AND VICE that are directly interpolated to the U and V grid. This should be fixed as these need to be rotated too
            if myvar in ['uice']:
                confM2R.grdMODEL.fieldSrc_u.data[:, :] = np.flipud(np.rot90(indata))
                field = confM2R.grdMODEL.regridSrc2Dst_u(confM2R.grdMODEL.fieldSrc_u, confM2R.grdMODEL.fieldDst_u)
            elif myvar in ['vice']:
                confM2R.grdMODEL.fieldSrc_v.data[:, :] = np.flipud(np.rot90(indata))
                field = confM2R.grdMODEL.regridSrc2Dst_v(confM2R.grdMODEL.fieldSrc_v, confM2R.grdMODEL.fieldDst_v)
            elif myvar in ['uvel']:
                confM2R.grdMODEL.fieldSrc_u.data[:, :] = np.flipud(np.rot90(indata))
                field = confM2R.grdMODEL.regridSrc2Dst_u(confM2R.grdMODEL.fieldSrc_u, confM2R.grdMODEL.fieldDst_rho)
            elif myvar in ['vvel']:
                confM2R.grdMODEL.fieldSrc_v.data[:, :] = np.flipud(np.rot90(indata))
                field = confM2R.grdMODEL.regridSrc2Dst_v(confM2R.grdMODEL.fieldSrc_v, confM2R.grdMODEL.fieldDst_rho)
            else:
                confM2R.grdMODEL.fieldSrc_rho.data[:, :] = np.flipud(np.rot90(indata))
                field = confM2R.grdMODEL.regridSrc2Dst_rho(confM2R.grdMODEL.fieldSrc_rho, confM2R.grdMODEL.fieldDst_rho)
            
            # Since ESMF uses coordinates (x,y) we need to rotate and flip to get back to (y,x) order.
            field = np.fliplr(np.rot90(field.data, 3))
           
        if confM2R.usefilter:
            field = laplacefilter(field, 1000, toxi, toeta)
            field = field * mymask
            
        array1[k, :, :] = field

    return array1

def setupIndexes(confM2R, myvar):
    if myvar in ["uice"]:
        indexROMS_Z_ST = (confM2R.grdMODEL.nlevels, confM2R.grdROMS.eta_u, confM2R.grdROMS.xi_u)
        toxi           = confM2R.grdROMS.xi_u
        toeta          = confM2R.grdROMS.eta_u
        mymask         = confM2R.grdROMS.mask_u
    elif myvar in ["vice"]:
        indexROMS_Z_ST = (confM2R.grdMODEL.nlevels, confM2R.grdROMS.eta_v, confM2R.grdROMS.xi_v)
        toxi           = confM2R.grdROMS.xi_v
        toeta          = confM2R.grdROMS.eta_v
        mymask         = confM2R.grdROMS.mask_v
    else:
        indexROMS_Z_ST = (confM2R.grdMODEL.nlevels, confM2R.grdROMS.eta_rho, confM2R.grdROMS.xi_rho)
        toxi           = confM2R.grdROMS.xi_rho
        toeta          = confM2R.grdROMS.eta_rho
        mymask         = confM2R.grdROMS.mask_rho
    return indexROMS_Z_ST, toxi, toeta, mymask


def setupESMFInterpolationWeights(confM2R):
    if confM2R.useesmf:
        print("regridSrc2Dst at RHO points")
        confM2R.grdMODEL.fieldSrc_rho      = ESMF.Field(confM2R.grdMODEL.esmfgrid, "fieldSrc", staggerloc=ESMF.StaggerLoc.CENTER)
        confM2R.grdMODEL.fieldDst_rho      = ESMF.Field(confM2R.grdROMS.esmfgrid, "fieldDst",staggerloc=ESMF.StaggerLoc.CENTER)
        confM2R.grdMODEL.regridSrc2Dst_rho = ESMF.Regrid(confM2R.grdMODEL.fieldSrc_rho, confM2R.grdMODEL.fieldDst_rho,regrid_method=ESMF.RegridMethod.BILINEAR,unmapped_action=ESMF.UnmappedAction.IGNORE)

        print("regridSrc2Dst at U points")
        confM2R.grdMODEL.fieldSrc_u        = ESMF.Field(confM2R.grdMODEL.esmfgrid_u, "fieldSrc", staggerloc=ESMF.StaggerLoc.CENTER)
        confM2R.grdMODEL.fieldDst_u        = ESMF.Field(confM2R.grdROMS.esmfgrid_u, "fieldDst_u",staggerloc=ESMF.StaggerLoc.CENTER)
        confM2R.grdMODEL.regridSrc2Dst_u   = ESMF.Regrid(confM2R.grdMODEL.fieldSrc_u, confM2R.grdMODEL.fieldDst_rho,regrid_method=ESMF.RegridMethod.BILINEAR,unmapped_action=ESMF.UnmappedAction.IGNORE)

        print("regridSrc2Dst at V points")
        confM2R.grdMODEL.fieldSrc_v        = ESMF.Field(confM2R.grdMODEL.esmfgrid_v, "fieldSrc", staggerloc=ESMF.StaggerLoc.CENTER)
        confM2R.grdMODEL.fieldDst_v        = ESMF.Field(confM2R.grdROMS.esmfgrid_v, "fieldDst_v",staggerloc=ESMF.StaggerLoc.CENTER)
        confM2R.grdMODEL.regridSrc2Dst_v   = ESMF.Regrid(confM2R.grdMODEL.fieldSrc_v, confM2R.grdMODEL.fieldDst_rho,regrid_method=ESMF.RegridMethod.BILINEAR,unmapped_action=ESMF.UnmappedAction.IGNORE)