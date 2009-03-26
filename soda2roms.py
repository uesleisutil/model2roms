import os, sys
from netCDF4 import Dataset
import numpy as np
import interp2D
from datetime import datetime, timedelta
from netCDF4 import num2date, date2num
import printObject
import interpolation as interp
import IOwrite
import plotData

#""" Get self made modules"""
#dir='/Users/trond/Projects/PyLIB'

#if os.path.isdir(dir):
#    sys.path.append(dir)
import date
import grd
import clim2bry
import barotropic
import IOinitial

__author__   = 'Trond Kristiansen'
__email__    = 'trond.kristiansen@imr.no'
__created__  = datetime(2008, 8, 15)
__modified__ = datetime(2008, 8, 19)
__modified__ = datetime(2009, 3, 11)
__version__  = "1.3"
__status__   = "Development"

def VerticalInterpolation(var,grdROMS,grdSODA):
    
    outINDEX_ST   = (grdROMS.Nlevels,grdROMS.eta_rho,grdROMS.xi_rho)
    outINDEX_U    = (grdROMS.Nlevels,grdROMS.eta_u,grdROMS.xi_u)
    outINDEX_UBAR = (grdROMS.eta_u,grdROMS.xi_u)
    outINDEX_V    = (grdROMS.Nlevels,grdROMS.eta_v,grdROMS.xi_v)
    outINDEX_VBAR = (grdROMS.eta_v,grdROMS.xi_v)
    
    if var=='temperature':
        data_ST=grdROMS.t[:,:,:]
    if var=='salinity':
        data_ST=grdROMS.s[:,:,:]
    if var=='vvel':
        data_U=grdROMS.u2[:,:,:]
        data_V=grdROMS.v2[:,:,:]
        
        
    if var=='salinity' or var=='temperature':
        print 'Interpolating vertically for %s with dimensions %s x %s'%(var,grdROMS.xi_rho,grdROMS.eta_rho)
        outdata=np.zeros((outINDEX_ST),dtype=np.float64, order='Fortran')
    
        outdata = interp.interpolation.dovertinter(np.asarray(outdata,order='Fortran'),
                                                       np.asarray(data_ST,order='Fortran'),
                                                       np.asarray(grdROMS.depth,order='Fortran'),
                                                       np.asarray(grdROMS.z_r,order='Fortran'),
                                                       np.asarray(grdSODA.z_r,order='Fortran'),
                                                       int(grdROMS.Nlevels),
                                                       int(grdSODA.Nlevels),
                                                       int(grdROMS.xi_rho),
                                                       int(grdROMS.eta_rho),
                                                       int(grdROMS.xi_rho),
                                                       int(grdROMS.eta_rho))
        
                                    
      
    if var=='vvel':
        print 'Interpolating vertically for %s with dimensions %s x %s'%(var,grdROMS.xi_u,grdROMS.eta_u)
        outdataU=np.zeros((outINDEX_U),dtype=np.float64,order='Fortran')
        outdataUBAR=np.zeros((outINDEX_UBAR),dtype=np.float64)

        outdataU = interp.interpolation.dovertinter(np.asarray(outdataU,order='Fortran'),
                                                       np.asarray(data_U,order='Fortran'),
                                                       np.asarray(grdROMS.depth,order='Fortran'),
                                                       np.asarray(grdROMS.z_r,order='Fortran'),
                                                       np.asarray(grdSODA.z_r,order='Fortran'),
                                                       int(grdROMS.Nlevels),
                                                       int(grdSODA.Nlevels),
                                                       int(grdROMS.xi_u),
                                                       int(grdROMS.eta_u),
                                                       int(grdROMS.xi_rho),
                                                       int(grdROMS.eta_rho))
   
        print 'Interpolating vertically for %s with dimensions %s x %s'%(var,grdROMS.xi_v,grdROMS.eta_v)
        outdataV=np.zeros((outINDEX_V),dtype=np.float64,order='Fortran')
        outdataVBAR=np.zeros((outINDEX_VBAR),dtype=np.float64)
         
        outdataV = interp.interpolation.dovertinter(np.asarray(outdataV,order='Fortran'),
                                                       np.asarray(data_V,order='Fortran'),
                                                       np.asarray(grdROMS.depth,order='Fortran'),
                                                       np.asarray(grdROMS.z_r,order='Fortran'),
                                                       np.asarray(grdSODA.z_r,order='Fortran'),
                                                       int(grdROMS.Nlevels),
                                                       int(grdSODA.Nlevels),
                                                       int(grdROMS.xi_v),
                                                       int(grdROMS.eta_v),
                                                       int(grdROMS.xi_rho),
                                                       int(grdROMS.eta_rho))
    
    if var=='temperature':
        grdROMS.t2[:,:,:]=outdata #*grdROMS.mask_rho
        #plotData.contourMap(grdROMS,grdSODA,np.squeeze(outdata[29,:,:]),"1",var)
    if var=='salinity':
        grdROMS.s2[:,:,:]=outdata #*grdROMS.mask_rho
        #plotData.contourMap(grdROMS,grdSODA,np.squeeze(outdata[29,:,:]),"1",var)
    if var=='vvel':
        grdROMS.u3[:,:,:]= outdataU #*grdROMS.mask_u
        
        
        outdataUBAR  = barotropic.velocity.ubar(np.asarray(outdataU,order='Fortran'),
                                                np.asarray(outdataUBAR,order='Fortran'),
                                                np.asarray(grdROMS.z_w,order='Fortran'),
                                                grdROMS.Nlevels,
                                                grdROMS.xi_u,
                                                grdROMS.eta_u,
                                                grdROMS.xi_rho,
                                                grdROMS.eta_rho)
        grdROMS.ubar = outdataUBAR
        
   
        grdROMS.v3[:,:,:]= outdataV #*grdROMS.mask_v
        outdataVBAR  = barotropic.velocity.vbar(np.asarray(outdataV,order='Fortran'),
                                                np.asarray(outdataVBAR,order='Fortran'),
                                                np.asarray(grdROMS.z_w,order='Fortran'),
                                                grdROMS.Nlevels,
                                                grdROMS.xi_v,
                                                grdROMS.eta_v,
                                                grdROMS.xi_rho,
                                                grdROMS.eta_rho)
        
        grdROMS.vbar = outdataVBAR
    
def HorizontalInterpolation(var,grdROMS,grdSODA,data):
    print 'Start %s horizontal interpolation for %s'%(grdSODA.grdType,var)
    
    if grdSODA.grdType=='regular':
        if var=='temperature':
            interp2D.doHorInterpolationRegularGrid(var,grdROMS,grdSODA,data)
        if var=='salinity':
            interp2D.doHorInterpolationRegularGrid(var,grdROMS,grdSODA,data)
        if var=='ssh':
            interp2D.doHorInterpolationSSHRegularGrid(var,grdROMS,grdSODA,data)
            grdROMS.ssh=grdROMS.ssh*grdROMS.mask_rho
        if var=='uvel':
            interp2D.doHorInterpolationRegularGrid(var,grdROMS,grdSODA,data)
        if var=='vvel':
            interp2D.doHorInterpolationRegularGrid(var,grdROMS,grdSODA,data)
            
    if grdSODA.grdType=='irregular':
        if var=='temperature':
            interp2D.doHorInterpolationIrregularGrid(var,grdROMS,grdSODA,temp)
        if var=='salinity':
            interp2D.doHorInterpolationIrregularGrid(var,grdROMS,grdSODA,salt)
        if var=='ssh':
            interp2D.doHorInterpolationSSHIrregularGrid(var,grdROMS,grdSODA,ssh)
        if var=='uvel':
            interp2D.doHorInterpolationIrregularGrid('uvel',grdROMS,grdSODA,uvel)
        if var=='vvel':
            interp2D.doHorInterpolationIrregularGrid('vvel',grdROMS,grdSODA,vvel)
    
    
    if var=='vvel':
        """
        First rotate the values of U, V at rho points with the angle, and then interpolate
        the rho point values to U and V points and save the result
        """
        
        urot=np.zeros((int(grdSODA.Nlevels),int(grdROMS.eta_rho),int(grdROMS.xi_rho)), np.float64)
        vrot=np.zeros((int(grdSODA.Nlevels),int(grdROMS.eta_rho),int(grdROMS.xi_rho)), np.float64)
       
        urot, vrot = interp.interpolation.rotate(np.asarray(urot,order='Fortran'),
                                                 np.asarray(vrot,order='Fortran'),
                                                 np.asarray(grdROMS.u,order='Fortran'),
                                                 np.asarray(grdROMS.v,order='Fortran'),
                                                 np.asarray(grdROMS.angle,order='Fortran'),
                                                 int(grdROMS.xi_rho),
                                                 int(grdROMS.eta_rho),
                                                 int(grdSODA.Nlevels))
        
        
        Zu=np.zeros((int(grdSODA.Nlevels),int(grdROMS.eta_u),int(grdROMS.xi_u)), np.float64)
        Zv=np.zeros((int(grdSODA.Nlevels),int(grdROMS.eta_v),int(grdROMS.xi_v)), np.float64)
        
        """
        Interpolate from RHO points to U and V points for velocities
        """
    
        Zu = interp.interpolation.rho2u(np.asarray(Zu,order='Fortran'),
                                        np.asarray(urot,order='Fortran'),
                                        int(grdROMS.xi_rho),
                                        int(grdROMS.eta_rho),
                                        int(grdSODA.Nlevels))
        
        grdROMS.u2[:,:,:]=Zu
        #plotData.contourMap(grdROMS,grdSODA,Zu[0,:,:],"1",'urot')
        
        Zv = interp.interpolation.rho2v(np.asarray(Zv,order='Fortran'),
                                        np.asarray(vrot,order='Fortran'),
                                        int(grdROMS.xi_rho),
                                        int(grdROMS.eta_rho),
                                        int(grdSODA.Nlevels))

        grdROMS.v2[:,:,:]=Zv
        #plotData.contourMap(grdROMS,grdSODA,Zv[0,:,:],"1",'vrot')

def getTime(grdROMS,grdSODA,year,ID):
    """
    Find the day and month that the SODA file respresents based on the year and ID number.
    Each SODA file represents a 5 day average, therefore we let the date we find be the first day
    of those 5 days. Thats the reason we subtract 4 below for day of month.
    """
    """
    Create a date object to keep track of Julian dates etc.
    Also create a reference date starting at 1948/01/01.
    Go here to check results:http://lena.gsfc.nasa.gov/lenaDEV/html/doy_conv.html
    """
    ref_date = date.Date()
    ref_date.day=1
    ref_date.month=1
    ref_date.year=1948
    jdref=ref_date.ToJDNumber()
    
    days=0.0; month=1;loop=True
    
    while loop is True:
        
        d=date.NumberDaysMonth(month,year)
        if days+d<int(ID)*5:
            days=days+d
            month+=1
        else:
            day=int(int(ID)*5-days)
            loop=False
            
    soda_date = date.Date()
    soda_date.day=day
    soda_date.month=month
    soda_date.year=year
    jdsoda=soda_date.ToJDNumber()
    
    grdROMS.time=(jdsoda-jdref)
    grdROMS.reftime=jdref
   
    print '\nCurrent time of SODA file : %s/%s/%s'%(soda_date.year,soda_date.month,soda_date.day)
    
def find_subset_indices(grdSODA,min_lat,max_lat,min_lon,max_lon):
    """
    Get the indices that covers the new grid, and enables us to only store a subset of
    the large input grid.
    """
    lat=grdSODA.lat[:,0]
    lon=grdSODA.lon[0,:]
    
    distances1 = []
    distances2 = []
    indices=[]
    index=0
    for point in lat:
        s1 = max_lat-point # (vector subtract)
        s2 = min_lat-point # (vector subtract)
        distances1.append((np.dot(s1, s1), point, index))
        distances2.append((np.dot(s2, s2), point, index))
        index=index+1
        
    distances1.sort()
    distances2.sort()
    indices.append(distances1[0])
    indices.append(distances2[0])
    
    distances1 = []
    distances2 = []
    index=0
   
    for point in lon:
        s1 = max_lon-point # (vector subtract)
        s2 = min_lon-point # (vector subtract)
        distances1.append((np.dot(s1, s1), point, index))
        distances2.append((np.dot(s2, s2), point, index))
        index=index+1
       
    distances1.sort()
    distances2.sort()
    indices.append(distances1[0])
    indices.append(distances2[0])
    
    """ Save final product: max_lat_indices,min_lat_indices,max_lon_indices,min_lon_indices"""
      
    grdSODA.minJ=indices[1][2]
    grdSODA.maxJ=indices[0][2]
    grdSODA.minI=indices[3][2]
    grdSODA.maxI=indices[2][2]
    
  
def convertSODA2ROMS(years,IDS,climName,initName,sodapath,romsgridpath):

    fileNameIn=sodapath+'SODA_2.0.2_'+str(years[0])+'_'+str(IDS[0])+'.cdf'
  
    """
    First time in loop, get the essential old grid information
    SODA data already at Z-levels. No need to interpolate to fixed depths,
    but we use the one we have
    """
    
    grdSODA = grd.grdClass(fileNameIn,"SODA")
    grdROMS = grd.grdClass(romsgridpath,"ROMS")
    
    """Now we want to subset the data to avoid storing more information than we need.
    We do this by finding the indices of maximum and minimum latitude and longitude in the matrixes"""
    
    find_subset_indices(grdSODA,min_lat=30, max_lat=90, min_lon=0, max_lon=360)
    #grdSODA.minJ=None
    #grdSODA.maxJ=None
    #grdSODA.minI=None
    #grdSODA.maxI=None
    grdSODA.lat=grdSODA.lat[grdSODA.minJ:grdSODA.maxJ,grdSODA.minI:grdSODA.maxI]
    grdSODA.lon=grdSODA.lon[grdSODA.minJ:grdSODA.maxJ,grdSODA.minI:grdSODA.maxI]
    
   
    print "\n---> Selected area in output file spans from (longitude=%3.2f,latitude=%3.2f) to (longitude=%3.2f,latitude=%3.2f)"%(grdROMS.lon_rho.min(),grdROMS.lat_rho.min(),grdROMS.lon_rho.max(),grdROMS.lat_rho.max())
    print "---> Selected area in input file spans from  (longitude=%3.2f,latitude=%3.2f) to (longitude=%3.2f,latitude=%3.2f)\n"%(grdSODA.lon.min(),grdSODA.lat.min(),grdSODA.lon.max(),grdSODA.lat.max())
    
    print '\n---> Finished initializing'
    print '\n--------------------------\n'
    
    time=0
    
    for year in years:
        
        firstRun = True ;
        
        for ID in IDS:
            file="SODA_2.0.2_"+str(year)+"_"+str(ID)+".cdf"
            filename=sodapath+file
    
            getTime(grdROMS,grdSODA,year,ID)
            
            cdf = Dataset(filename)
            
            """Each SODA file consist only of one time step. Get the subset data selected, and
            store that time step in a new array:"""
            
        
            if firstRun is True:
                firstRun = False
                indexTMP_ST    = (grdSODA.Nlevels,grdROMS.eta_rho,grdROMS.xi_rho)
                indexTMP_U     = (grdSODA.Nlevels,grdROMS.eta_u,grdROMS.xi_u)
                indexTMP_V     = (grdSODA.Nlevels,grdROMS.eta_v,grdROMS.xi_v)
                
                indexSODA_Z    = (grdSODA.Nlevels,len(grdSODA.lat),len(grdSODA.lon))
                indexROMS_Z_ST = (grdSODA.Nlevels,grdROMS.eta_rho,grdROMS.xi_rho)
                indexROMS_S_ST = (grdROMS.Nlevels,grdROMS.eta_rho,grdROMS.xi_rho)
                indexROMS_SSH  = (grdROMS.eta_rho,grdROMS.xi_rho)
                
                indexROMS_Z_U = (grdSODA.Nlevels,grdROMS.eta_u,grdROMS.xi_u)
                indexROMS_S_U = (grdROMS.Nlevels,grdROMS.eta_u,grdROMS.xi_u)
                
                indexROMS_Z_V = (grdSODA.Nlevels,grdROMS.eta_v,grdROMS.xi_v)
                indexROMS_S_V = (grdROMS.Nlevels,grdROMS.eta_v,grdROMS.xi_v)
                indexROMS_UBAR = (grdROMS.eta_u,grdROMS.xi_u)
                indexROMS_VBAR = (grdROMS.eta_v,grdROMS.xi_v)
                
        
            """
            All variables for all time are now stored in arrays. Now, start the interpolation to the
            new grid for all variables and then finally write results to file.
            """
            vars=['temperature','salinity','ssh','uvel','vvel']
            #vars=['temperature']
            #vars=['uvel','vvel']
            #vars=['ssh']
            for var in vars:
                if var=='temperature':
                    data = np.array(cdf.variables["TEMP"][0,:,grdSODA.minJ:grdSODA.maxJ,grdSODA.minI:grdSODA.maxI])
                    grdROMS.t=np.zeros((indexROMS_Z_ST),dtype=np.float64)
                    data_ST =np.zeros((indexTMP_ST),dtype=np.float64)
                    grdROMS.t2=np.zeros((indexROMS_S_ST),dtype=np.float64)
                   
                    
                if var=='salinity':
                    data = np.array(cdf.variables["SALT"][0,:,grdSODA.minJ:grdSODA.maxJ,grdSODA.minI:grdSODA.maxI])
                    grdROMS.s=np.zeros((indexROMS_Z_ST),dtype=np.float64)
                    data_ST =np.zeros((indexTMP_ST),dtype=np.float64)
                    grdROMS.s2=np.zeros((indexROMS_S_ST),dtype=np.float64)
                  
                    
                if var=='ssh':
                    data = np.array(cdf.variables["SSH"][0,grdSODA.minJ:grdSODA.maxJ,grdSODA.minI:grdSODA.maxI])
                    grdROMS.ssh=np.zeros((indexROMS_SSH),dtype=np.float64)
                   
                    
                if var=='uvel':
                    data = np.array(cdf.variables["U"][0,:,grdSODA.minJ:grdSODA.maxJ,grdSODA.minI:grdSODA.maxI])
                    grdROMS.u=np.zeros((indexROMS_Z_ST),dtype=np.float64)
                    data_U  =np.zeros((indexTMP_U),dtype=np.float64)
                    grdROMS.u2=np.zeros((indexROMS_Z_U),dtype=np.float64)
                    grdROMS.u3=np.zeros((indexROMS_S_U),dtype=np.float64)
                   
                    
                if var=='vvel':
                    data     = np.array(cdf.variables["V"][0,:,grdSODA.minJ:grdSODA.maxJ,grdSODA.minI:grdSODA.maxI])
                    grdROMS.v=np.zeros((indexROMS_Z_ST),dtype=np.float64)
                    data_V   =np.zeros((indexTMP_V),dtype=np.float64)
                    grdROMS.v2=np.zeros((indexROMS_Z_V),dtype=np.float64)
                    grdROMS.v3=np.zeros((indexROMS_S_V),dtype=np.float64)
                    
                     
                HorizontalInterpolation(var,grdROMS,grdSODA,data)
                VerticalInterpolation(var,grdROMS,grdSODA)
                
                
                IOwrite.writeClimFile(grdROMS,time,climName,var)
            
            if time==0:
                IOinitial.createInitFile(grdROMS,time,initName,var)
                
            cdf.close()    
            time+=1
          
            
          
            
    