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

dir = './from_runs/'
prefix='china_boff_hi'
mns=['20141201', '20150101', '20150201']
#mns=['1','2','3','4','5','6','7','8','9','10','11','12']
nbins = 15 # number of bins
ibins = nbins
ndry = 7 # number of dry aerosol components
molwgt = [96., 58.5, 12., 12., 12., 12.,100.,18.] # mol wgt [g/mol]
#xk = [1.e-21*2.**(-10)] # bin lower limits [kg]
#for k in range(0,nbins):

## GEOS-3 reduced grid
#nlevs=30
#etae=np.array([0.000000,0.000056,0.000201,0.000607,0.001641,0.004067,\
#      0.009283,0.019782,0.028358,0.040166,0.056379,0.078503,\
#      0.092357,0.108655,0.127828,0.150385,0.176922,0.208144,\
#      0.244867,0.288076,0.338908,0.398712,0.469063,0.551839,\
#      0.649201,0.744379,0.831017,0.903295,0.955994,0.985110,\
#      1.000000]) # box edge sigma coordinates
#etam=np.array([0.0000280,0.0001290,0.0004040,0.0011240,0.0028540,0.0066750,\
#      0.0145320,0.0240700,0.0342620,0.0482720,0.0674410,0.0854300,\
#      0.1005060,0.1182410,0.1391060,0.1636530,0.1925330,0.2265060,\
#      0.2664720,0.3134920,0.3688100,0.4338870,0.5104510,0.6005200,\
#      0.6967900,0.7876980,0.8671560,0.9296450,0.9705520,0.9925550])# box middle
#psurf=1000 #mb
#ptop=0.01 #mb

# GEOS-5 reduced grid
nlevs=47
etae=np.array([0.000000,0.000055,0.000199,0.000601,0.001625,0.004026,\
      0.009191,0.019586,0.028077,0.039768,0.055820,0.077726,\
      0.091442,0.107578,0.126563,0.148896,0.175170,0.206167,\
      0.242774,0.285974,0.335486,0.373114,0.410759,0.448431,\
      0.486118,0.523819,0.561527,0.599251,0.636974,0.674708,\
      0.699867,0.725026,0.750186,0.775350,0.800515,0.820648,\
      0.835748,0.850848,0.865949,0.881051,0.896152,0.911253,\
      0.926356,0.941459,0.956562,0.971665,0.986769,1.001796])# box edge sigma coordinates
etam=np.array([0.000028,0.000127,0.000400,0.001113,0.002825,0.006609,\
      0.014389,0.023832,0.033923,0.047794,0.066773,0.084584,\
      0.099510,0.117070,0.137729,0.162033,0.190668,0.224471,\
      0.264374,0.310730,0.354300,0.391937,0.429595,0.467274,\
      0.504968,0.542673,0.580389,0.618113,0.655841,0.687287,\
      0.712447,0.737606,0.762768,0.787933,0.810582,0.828198,\
      0.843298,0.858399,0.873500,0.888601,0.903703,0.918805,\
      0.933908,0.949010,0.964113,0.979217,0.994283]) # box middle
psurf=1000 #mb
ptop=0.01 #mb

sdhead = ['NK','SF','SS','ECIL','ECOB','OCIL','OCOB','DUST','AW']
#sdhead = ['NK']

###################
# Main code
###################

CN3=[]
CN10=[]
CN40=[]
CN80=[]
CN100=[]
H2SO4=[]
OH=[]
nuc=[]
OC=[]
BC=[]

# make arrays of pressure
prese=[]
for ll in range(0,nlevs+1):
   l=nlevs-ll-1
   prese.append(ptop+etae[l]*(psurf-ptop))
prese=np.array(prese)

presm=[]
for ll in range(0,nlevs):
   l=nlevs-ll-1
   presm.append(ptop+etam[l]*(psurf-ptop))
presm=np.array(presm)
data=[]
dist=[]
for d in range(0,len(mns)):
   fname=dir+prefix+'_'+mns[d]+'.nc'
   #fnameOH=dir+prefix+'_'+mns[d]+'_CHEM-L.nc'
   #fnameNuc=dir+prefix+'_'+mns[d]+'_nuc.nc'
   #fname=dir+prefix+'.nc'
   nc_file = Dataset(fname,'r')
   #nc_fileOH = Dataset(fnameOH,'r')
   #nc_fileNuc = Dataset(fnameNuc,'r')
   data = []
   for c in range(0,len(sdhead)):
      sizedist = []
      for k in range(0,ibins):
         sizedist.append(nc_file.variables['IJ_AVG_S__'+sdhead[c]+str(k+1)][:])
      data.append(sizedist)
   dist.append(data)
   #data=np.array(data)
   #H2SO4.append(nc_file.variables['IJ_AVG_S__H2SO4'][0,:,:,:])
   #OH.append(nc_fileOH.variables['CHEM_L_S__OH'][0,:,:,:])
   #nuc.append(nc_fileNuc.variables['TOMAS_3D__NUC_1NM'][0,:,:,:])
   lat=nc_file.variables['lat'][:]
   lon=nc_file.variables['lon'][:]
   nc_file.close()



dist = np.array(dist)

lat = np.array(lat)
lon = np.array(lon)
dist = np.array(dist)[:,:,:,0,:,:,:] #AH: added [:,:,0,:,:,:] as JK did that to get rid of the unused 3rd dimension 
np.savez(prefix+'_dist_w1415.npz', dist=dist, la=lat,lo=lon)


