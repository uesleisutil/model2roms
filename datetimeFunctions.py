from netCDF4 import Dataset, date2num, num2date
import os
from datetime import datetime, timedelta
import calendar

"""
Created by Trond Kristiansen
https://github.com/trondkr/model2roms
"""

def createlistofmonths(confM2R, currentyear):

    if currentyear == confM2R.startdate.year:
        IDS = [confM2R.startdate.month + m for m in range(13 - confM2R.startdate.month)]
        # months from start month to end of first year

    elif currentyear == confM2R.enddate.year:
        IDS = [1 + m for m in range(int(confM2R.enddate.month))]
        # months from first month to last month of last year

    elif confM2R.startdate.year < currentyear < confM2R.enddate.year:
        # months from first month to last month of last year
        IDS = [1 + m for m in range(12)]

    if confM2R.startdate.year == confM2R.enddate.year:
        # months from first month to last month of last year
        IDS = [confM2R.startdate.month + m for m in range(confM2R.enddate.month - confM2R.startdate.month + 1)]

    if not IDS:
        print("Unable to generate IDS for time looping: main.py -> func createIDS")
        sys.exit()

    return IDS


def createlistofdays(confM2R, year, month):
    days = []
    if confM2R.timefrequencyofinputdata == 'day':
        daystep = 5
        if daystep > 1:
            print("WARNING!")
            print("----------------------------------------------------------------------")
            print("You are only using every {} days of input data! (model2roms.py)".format(daystep))
            print("----------------------------------------------------------------------")

        # Regulary we want all days in each month                
        ndays = int(calendar.monthrange(year, month)[1])
        days = [d + 1 for d in range(0, ndays, daystep)]

        # Exceptions:
        # We start in the first month after day one
        if (month == confM2R.startdate.month and year == confM2R.startdate.year and confM2R.startdate.day > 1):
            days = [i for i in range(confM2R.startdate.day, ndays, daystep)]

        # We finish in the last month before last day of month
        if (month == confM2R.enddate.month and year == confM2R.enddate.year):
            days = [i + 1 for i in range(0, confM2R.enddate.day, daystep)]

        # We start and end on different days in the same month         
        if (confM2R.startdate.month == confM2R.enddate.month and confM2R.startdate.year == confM2R.enddate.year):
            days = [i for i in range(confM2R.startdate.day, confM2R.enddate.day, daystep)]
        
        return days

    if confM2R.timefrequencyofinputdata == 'month':
        days = [15]

    return days
