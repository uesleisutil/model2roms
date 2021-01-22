import time
from datetime import datetime, timedelta
import configM2R
import model2roms
import clim2bry
import atmosForcing

"""
Created by Trond Kristiansen
https://github.com/trondkr/model2roms
"""

def run():
    confM2R = configM2R.Model2romsConfig()

    if confM2R.createatmosforcing or confM2R.createoceanforcing:

        if confM2R.createoceanforcing:
            model2roms.convertMODEL2ROMS(confM2R)

            clim2bry.writebry(confM2R)

      #  if confM2R.createAtmosForcing:
      #      atmosForcing.createAtmosFileUV(confM2R)

    print('Finished ' + time.ctime(time.time()))

run()
