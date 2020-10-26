# This script calculates the CN and CCN global maps from nc files

###################
# Modules
###################

import datetime as dt
import numpy as np
#from scipy.io import netcdf
from netCDF4 import Dataset
import sys
###################
# Inputs
###################

dir = '/pierce-scratch/kbilsback/geos_chem/v12run_dirs_shipping/geosfp_4x5_TOMAS15_noscale/OutputDir/'
prefix='shipping'
newfile = 'ship_noscale_aw'
mns=['20130101', '20130201', '20130301', '20130401', '20130501', '20130601',
     '20130701', '20130801', '20130901', '20131001', '20131101', '20131201']

ibins = 15

sdhead = ['NK','SF','SS','ECIL','ECOB','OCIL','OCOB','DUST','AW']

###################
# Main code
###################

data=[]
dist= np.zeros((len(mns), len(sdhead) + 2, ibins, 47, 46, 72))
nit=[]
nh4=[]

for d in range(0,len(mns)):
   fname=dir+prefix+'_'+mns[d]+'.nc'
   nc_file = Dataset(fname,'r')
   nit.append(nc_file.variables['IJ_AVG_S__NIT'][0])
   nh4.append(nc_file.variables['IJ_AVG_S__NH4'][0])
   data = []
   for c in range(0,len(sdhead)):
      sizedist = []
      for k in range(0,ibins):
         sizedist.append(nc_file.variables['IJ_AVG_S__'+sdhead[c]+str(k+1)][0])
      data.append(sizedist)
   dist[d,:9] = data
   lat=nc_file.variables['lat'][:]
   lon=nc_file.variables['lon'][:]
   nc_file.close()

dist = np.array(dist)

# distribution ammonia and nitrate to sulfate
aw = dist[:,8]
aw_total = np.sum(aw,1)

nit_dist = np.zeros((len(mns),ibins, 47, len(lat), len(lon)))
nh4_dist = np.zeros((len(mns),ibins, 47, len(lat), len(lon)))
for i in range(ibins):
   nit_dist[:,i] = nit * (aw[:,i] / aw_total)
   nh4_dist[:,i] = nh4 * (aw[:,i] / aw_total)

# construct final array
lat = np.array(lat)
lon = np.array(lon)
dist[:,9] = nit_dist
dist[:,10] = nh4_dist

np.savez(newfile + '_dist.npz', dist=dist, la=lat,lo=lon)


