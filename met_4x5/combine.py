# Combine met file and albedo files for AIE calculations
# Kelsey Bilsback
# January 2019

##################################
# Import libraries
##################################

from pylab import *
from netCDF4 import Dataset

###################################################
# Read existing npz files with met and albedo data
###################################################

# load cloud fraction
f1 = load('clfrac_4x5_w1415.npz')
clfrac = f1['clfrac']
f1.close()

# load cloud ice water mixing ratios
f2 = load('qi_4x5_w1415.npz')
qi = f2['qi']
f2.close()

# load cloud liquid water mixing ratios
f3 = load('ql_4x5_w1415.npz')
ql = f3['ql']
f3.close()

# load specific humidity
f4 = load('qv_4x5_w1415.npz')
qv = f4['qv']
f4.close()

# load surface pressure
f5 = load('surf_press_4x5_w1415.npz')
surf_press = f5['surf_press']
f5.close()

# load temperature
f6 = load('temperature_4x5_w1415.npz')
temperature = f6['temperature']
f6.close()

# load direct and diffuse surface albedos at vis and NIR wavelengths
f7 = load('seaice_albedo_4x5.npz')
albdrvis = f7['albdrvis'][:]
albdfvis = f7['albdfvis'][:]
albdrnir = f7['albdrnir'][:]
albdfnir = f7['albdfnir'][:]
lat = f7['lat'][:]
lon = f7['lon'][:]
f7.close()

##################################
# Write netCDF file
##################################

# write grid file
nc_w_fid = Dataset('grid_4x5.nc', 'w', format='NETCDF4')
nc_w_fid.description = 'Lat and lon grid'
nc_w_fid.history = 'Created ' + datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

nc_w_fid.createDimension('lat', 46)
nc_w_fid.createDimension('lon', 72)

lat_w = nc_w_fid.createVariable('lat', np.float32, ('lat',))
lon_w = nc_w_fid.createVariable('lon', np.float32, ('lon',))

lat_w[:] = lat
lon_w[:] = lon

nc_w_fid.close()

# file descriptors
nc_w_fid = Dataset('met_albedo_4x5_w1415.nc', 'w', format='NETCDF4')
nc_w_fid.description = 'Met and albedo data for AIE calculation'
nc_w_fid.history = 'Created ' + datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

# define file dimensions
nc_w_fid.createDimension('time', None)
nc_w_fid.createDimension('lat', 46)
nc_w_fid.createDimension('lon', 72)
nc_w_fid.createDimension('lev', 47)

# create identify variables
time_w = nc_w_fid.createVariable('time', np.float32, ('time',))
lev_w = nc_w_fid.createVariable('lev', np.float32, ('lev',))
lat_w = nc_w_fid.createVariable('lat', np.float32, ('lat',))
lon_w = nc_w_fid.createVariable('lon', np.float32, ('lon',))

# create albedo and met variables
albdrvis_w = nc_w_fid.createVariable('albdrvis', np.float32, ('time','lat','lon',))
albdfvis_w = nc_w_fid.createVariable('albdfvis', np.float32, ('time','lat','lon',))
albdrnir_w = nc_w_fid.createVariable('albdrnir', np.float32, ('time','lat','lon',))
albdfnir_w = nc_w_fid.createVariable('albdfnir', np.float32, ('time','lat','lon',))

surf_press_w = nc_w_fid.createVariable('surf_press', np.float32, ('time','lat','lon',))
ql_w = nc_w_fid.createVariable('ql', np.float32, ('time','lev','lat','lon',))
qi_w = nc_w_fid.createVariable('qi', np.float32, ('time','lev','lat','lon',))
qv_w = nc_w_fid.createVariable('qv', np.float32, ('time','lev','lat','lon',))
temperature_w = nc_w_fid.createVariable('temperature', np.float32, ('time','lev','lat','lon',))
clfrac_w = nc_w_fid.createVariable('clfrac', np.float32, ('time','lev','lat','lon',))

# fill in identity variables
lev_w[:] = arange(0, 47, 1)
lat_w[:] = lat
lon_w[:] = lon
time_w[:] = arange(0, 12, 1)

# fill in albedo variables
albdrvis_w[:,:,:] = albdrvis
albdfvis_w[:,:,:] = albdfvis
albdrnir_w[:,:,:] = albdrnir
albdfnir_w[:,:,:] = albdfnir

# fill in met variables
surf_press_w[:,:,:] = surf_press
ql_w[:,:,:,:] = ql[:,:47,:,:]
qi_w[:,:,:,:] = qi[:,:47,:,:]
qv_w[:,:,:,:] = qv[:,:47,:,:]
temperature_w[:,:,:,:] = temperature[:,:47,:,:]
clfrac_w[:,:,:,:] = clfrac[:,:47,:,:]

# add variable attributes
time_w.setncatts({'units':'months since 2014-12',
                  'long_name':'Time'})

lat_w.setncatts({'units':'degrees_north',
                 'long_name':'Latitude'})

lon_w.setncatts({'units':'degrees_east',
                 'long_name':'Longitude'})

lev_w.setncatts({'units':'unitless',
                 'long_name':'Levels'})

albdrvis_w.setncatts({'units':'unitless',
                      'long_name':'Direct surface albedo (visible)'})

albdrnir_w.setncatts({'units':'unitless',
                      'long_name':'Direct surface albedo (near-IR)'})

albdfvis_w.setncatts({'units':'unitless',
                      'long_name':'Diffuse  surface albedo (visible)'})

albdfnir_w.setncatts({'units':'unitless',
                      'long_name':'Diffuse surface albedo (near-IR)'})

surf_press_w.setncatts({'units':'hPa',
                        'long_name':'Surface pressure'})

temperature_w.setncatts({'units':'K',
                         'long_name':'Temperature'})

qv_w.setncatts({'units':'kg kg-1',
                'long_name':'Specific humidity'})

ql_w.setncatts({'units':'kg kg-1',
                'long_name':'Cloud liquid water mixing ratio'})

qi_w.setncatts({'units':'kg kg-1',
                'long_name':'Cloud ice water mixing ratio'})
                
clfrac_w.setncatts({'units':'unitless',
                    'long_name':'Total cloud fraction in grid box'})

nc_w_fid.close()
