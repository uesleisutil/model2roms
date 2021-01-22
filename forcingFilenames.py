import os

# Main function called from model2roms
def getFilename(confM2R,year,month,defaultvar):
    if confM2R.indatatype == 'SODA3':
        if defaultvar is None:defaultvar="salinity"
        filenamein = getSODA3filename(confM2R, year, month, defaultvar)
    if confM2R.indatatype == 'GLORYS':
        if defaultvar is None:defaultvar="temperature"
        filenamein = getGLORYSfilename(confM2R, year, month, defaultvar)
    return filenamein
    
# private functions called from within module
def getGLORYSfilename(confM2R, year, month, myvar):
    filename = confM2R.modelpath + 'glorys.nc'

    return filename

#def getSODA3filename(confM2R, year, month, myvar):
    if (myvar in ['cn', 'hi', 'hs']):
        return confM2R.modelpath + "soda.nc"
    else:
        return confM2R.modelpath + "soda.nc"