from pylab import *
import numpy as np
from netCDF4 import Dataset
import sys
import os

os.chdir('/shared-scratch/GEOS-Chem/ExtData/GEOS_4x5/GEOS_FP/')

start=[20160801, 20160901, 20161001, 20161101, 20161201, 20170101,
              20170201, 20170301, 20170401, 20170501, 20170601, 20170701]

dstart = [1, 1, 1, 1, 1, 1,
                    1, 1, 1, 1, 1, 1]
days =   [31, 30, 31, 30, 31, 31,
                    28, 31, 30, 31, 30, 31]

monthly_data = []
for mn in range(0, len(start)): #iterate through 12 months
    daily_data = []
    year = str(start[mn])[:4]
    month = str(start[mn])[4:6]
    os.chdir(year + '/' + month)
    for dy in range(0, days[mn]): #iterate through each day of the month
        date = start[mn] + dy # get date
        fname = 'GEOSFP.' + str(date)+'.A3cld.4x5.nc'
        print(fname)
        nc_file = Dataset(fname)
        nc_data = nc_file.variables['QL'][:]
        daily_data.append(np.mean(nc_data, axis=0))
    monthly_data.append(np.array(daily_data).mean(0))
    os.chdir('/shared-scratch/GEOS-Chem/ExtData/GEOS_4x5/GEOS_FP/')
    
monthly_data = array(monthly_data)

os.chdir('/pierce-scratch/kbilsback/geos_chem/projects/atom/radiative_effects_4x5/met_4x5')

np.savez('ql_4x5_20160801.npz', ql=monthly_data)
