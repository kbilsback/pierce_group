# Modify albedo to include sea ice

from pylab import *
from netCDF4 import Dataset
#from tomas_class import TOMAS as tom

#### -------------------------------------------------------------------------------####
#### Assign albedos
#### -------------------------------------------------------------------------------####
alpha_ocean = 0.07
alpha_ice_vis = 0.73
alpha_ice_nir = 0.33
alpha_snow_vis = 0.96
alpha_snow_nir = 0.68
#http://www.cesm.ucar.edu/models/atm-cam/docs/description/node35.html

#### -------------------------------------------------------------------------------####
#### Read in sea ice fraction from MERRA
#### -------------------------------------------------------------------------------####
f1 = load('sea_ice_fraction_climatology.npz')
FRSICE = f1['FRSICE']
FRSICE = np.concatenate((FRSICE[7:12], FRSICE[0:7]),0)
FRSNOW = f1['FRSNOW']
FRSNOW = np.concatenate((FRSNOW[7:12], FRSNOW[0:7]),0)
lat = f1['lat']
lon = f1['lon']
f1.close()

#### -------------------------------------------------------------------------------####
#### Read in MODIS albedo file
#### -------------------------------------------------------------------------------####
f1 = load('modis_albedo_4x5.npz')

# ALBEDO
ALBDRVIS = f1['albdrvis']
ALBDRVIS = np.concatenate((ALBDRVIS[7:12], ALBDRVIS[0:7]),0)
ALBDFVIS = f1['albdfvis']
ALBDFVIS = np.concatenate((ALBDFVIS[7:12], ALBDFVIS[0:7]),0)
ALBDRNIR = f1['albdrnir']
ALBDRNIR = np.concatenate((ALBDRNIR[7:12], ALBDRNIR[0:7]),0)
ALBDFNIR = f1['albdfnir']
ALBDFNIR = np.concatenate((ALBDFNIR[7:12], ALBDFNIR[0:7]),0)

# Crop to April -- should re-do all
#ALBDRVIS = ALBDRVIS[:,:,:]
#ALBDFVIS = ALBDFVIS[:,:,:]
#ALBDRNIR = ALBDRNIR[3,:,:]
#ALBDFNIR = ALBDFNIR[3,:,:]

#### -------------------------------------------------------------------------------####
#### Calculate new albedo
#### -------------------------------------------------------------------------------####
#FRSICE = ma.masked_less(FRSICE, 0.85).filled(0)

#alphaSnowIce_vis = FRSNOW*alpha_snow_vis + (1-FRSNOW)*alpha_ice_vis

ALBDRVIS_mod = np.where(ALBDRVIS <= 0.1, FRSICE*alpha_ice_vis + (1-FRSICE)*ALBDRVIS, ALBDRVIS)
ALBDFVIS_mod = np.where(ALBDFVIS <= 0.1, FRSICE*alpha_ice_vis + (1-FRSICE)*ALBDFVIS, ALBDFVIS)

ALBDRNIR_mod = np.where(ALBDRNIR <= 0.1, FRSICE*alpha_ice_nir + (1-FRSICE)*ALBDRNIR, ALBDRNIR)
ALBDFNIR_mod = np.where(ALBDFNIR <= 0.1, FRSICE*alpha_ice_nir + (1-FRSICE)*ALBDRNIR, ALBDFNIR)

#### -------------------------------------------------------------------------------####
#### Save values
#### -------------------------------------------------------------------------------####
np.savez('seaice_albedo_4x5.npz', ALBDRVIS_mod=ALBDRVIS_mod, ALBDFVIS_mod=ALBDFVIS_mod\
         ,ALBDRNIR_mod=ALBDRNIR_mod, ALBDFNIR_mod=ALBDFNIR_mod, lat=lat, lon=lon)

#### -------------------------------------------------------------------------------####
#### TEST PLOT
#### -------------------------------------------------------------------------------####
#close('all')

### Mod albedo
#levels = [0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1.0]
#figure()
#m,x,y,datap = tom.make_map_arctic(lat, lon, ALBDRVIS_mod[3,:,:])
#overlay = contourf(x,y,datap, levels)
#title('Modified Albedo Vis Jan')
#colorbar()
#savefig('alb_mod_vis.png')
#show()

#figure()
#m,x,y,datap = tom.make_map_arctic(lat, lon, ALBDRVIS[3,:,:])
#overlay = contourf(x,y,datap, levels)
#title('Original Albedo Vis Jan')
#colorbar()
#savefig('alb_orig_vis.png')
#show()

#figure()
#m,x,y,datap = tom.make_map_arctic(lat, lon, ALBDFVIS_mod[7,:,:])
#overlay = contourf(x,y,datap, levels)
#title('Modified Albedo Vis Aug')
#colorbar()
#savefig('alb_mod_nir.png')
#show()

#figure()
#m,x,y,datap = tom.make_map_arctic(lat, lon, ALBDFVIS[7,:,:])
#overlay = contourf(x,y,datap, levels)
#title('Original Albedo Vis Aug')
#colorbar()
#savefig('alb_org_nir.png')
#show()
