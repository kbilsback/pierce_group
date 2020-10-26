# This script makes an input file for RRTMG AIE

from pylab import *
from netCDF4 import Dataset
import sys

#### -------------------------------------------------------------------------------####
#### Routines
#### -------------------------------------------------------------------------------####
# This just returns monthly solar declination angle
def getMonSolDecAng():
    # Get Monthly Avg Power
    # solar declination angle
    daymonth=np.array([31.,28.,31.,30.,31.,30.,31.,31.,30.,31.,30.,31.])
    monSolDecAng=[] # monthly solar declination angle (radians)
    tot=0.
    for m in range(0,12):
        avgday=tot+daymonth[m]/2.
        monSolDecAng.append(0.409*np.cos(2.*np.pi*(avgday-173.)/365.))
        tot=tot+daymonth[m]
    monSolDecAng = np.array(monSolDecAng)
    return monSolDecAng

# ap [hPa] for 47 levels (48 edges)
ap = np.array([0.000000e+00, 4.804826e-02, 6.593752e+00, 1.313480e+01,
               1.961311e+01, 2.609201e+01, 3.257081e+01, 3.898201e+01,
               4.533901e+01, 5.169611e+01, 5.805321e+01, 6.436264e+01,
               7.062198e+01, 7.883422e+01, 8.909992e+01, 9.936521e+01,
               1.091817e+02, 1.189586e+02, 1.286959e+02, 1.429100e+02,
               1.562600e+02, 1.696090e+02, 1.816190e+02, 1.930970e+02,
               2.032590e+02, 2.121500e+02, 2.187760e+02, 2.238980e+02,
               2.243630e+02, 2.168650e+02, 2.011920e+02, 1.769300e+02,
               1.503930e+02, 1.278370e+02, 1.086630e+02, 9.236572e+01,
               7.851231e+01, 5.638791e+01, 4.017541e+01, 2.836781e+01,
               1.979160e+01, 9.292942e+00, 4.076571e+00, 1.650790e+00,
               6.167791e-01, 2.113490e-01, 6.600001e-02, 1.000000e-02])

# bp [unitless] for 47 levels (48 edges)
bp = np.array([1.000000e+00, 9.849520e-01, 9.634060e-01, 9.418650e-01,
               9.203870e-01, 8.989080e-01, 8.774290e-01, 8.560180e-01,
               8.346609e-01, 8.133039e-01, 7.919469e-01, 7.706375e-01,
               7.493782e-01, 7.211660e-01, 6.858999e-01, 6.506349e-01,
               6.158184e-01, 5.810415e-01, 5.463042e-01, 4.945902e-01,
               4.437402e-01, 3.928911e-01, 3.433811e-01, 2.944031e-01,
               2.467411e-01, 2.003501e-01, 1.562241e-01, 1.136021e-01,
               6.372006e-02, 2.801004e-02, 6.960025e-03, 8.175413e-09,
               0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
               0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
               0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
               0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00])

#### -------------------------------------------------------------------------------####
#### File names and options
#### -------------------------------------------------------------------------------####
# CDNC FILES -- made from CDNC.py
fname1 = 'CDNC_output/aw_ship_so2ocbcnoxscale_cdnc.npz'
fname2 = 'CDNC_output/aw_ship_noscale_cdnc.npz'

# OUTFILE NAME
fout ='aw_ship_so2ocbcnoxscale.nc'
CONTROL=False

# Only do rad-transfer on first 30 levels
ltop = 30
#### -------------------------------------------------------------------------------####
#### Get MET variables
#### -------------------------------------------------------------------------------####
# READ IN MET FILE
nc_file = Dataset('../met_4x5/met_albedo_4x5.nc')
lat = nc_file.variables['lat'][:]
lon = nc_file.variables['lon'][:]

LLPAR = 47
nmonths = 12

# Pressure Endge
surf_press = nc_file.variables['surf_press'][:]  #[time,lat,lon]
# Create 3d pressure array
PLEV = np.zeros((nmonths,LLPAR+1,len(lat),len(lon)))
for l in range(0, LLPAR+1):
    PLEV[:,l,:,:] = ap[l] + ( bp[l] * surf_press[:,:,:])

PCENTER = np.zeros((nmonths,LLPAR,len(lat),len(lon)))
for l in range(0, LLPAR):
    PCENTER[:,l,:,:] = (PLEV[:,l] + PLEV[:,l+1]) / 2

# TEMPERATURE
TLAY = nc_file.variables['temperature'][:]
TSFC = TLAY[:,0,:,:]

TLEV = np.zeros((12,LLPAR+1,len(lat),len(lon)))  #FIX THIS IN THE FUTURE
TLEV[:,:-1,:,:] = TLAY[:,:,:,:]
TLEV[:,-1,:,:] = TLAY[:,-1,:,:]
TLAY = (TLEV[:,:-1,:,:]+TLEV[:,1:,:,:])/2.

# ALBEDO
ALBDRVIS = nc_file.variables['albdrvis'][:]
ALBDFVIS = nc_file.variables['albdfvis'][:]
ALBDRNIR = nc_file.variables['albdrnir'][:]
ALBDFNIR = nc_file.variables['albdfnir'][:]

# CLOUD FRACTION
clfrac = nc_file.variables['clfrac'][:]

# SOLAR DECLINATION ANGLE - LOOK INTO THIS
dhrs= 3
HRS = 24//dhrs
monSolDecAng = getMonSolDecAng()
RA = np.array([19.42, 21.42, 23.39, 1.32, 3.15, 5.58, 7.28, 9.41, 11.34, 13.27, 15.16, 17.4])
LONS, LATS = np.meshgrid(lon, lat)
SUNCOS = []
SUNCOSHA = []

for i in range(0, 12):
    for j in range(0,HRS):
        k = j*dhrs*360/24.
        HA = k - LONS - RA[i]
        
        SUNCOSHA.append(sin(LATS*pi/180.)*sin(monSolDecAng[i]) + cos(LATS*pi/180.)*cos(monSolDecAng[i])*cos(HA*pi/180.))
    SUNCOS.append(SUNCOSHA)
    SUNCOSHA = []

SUNCOS = ma.masked_less(array(SUNCOS), 0.).filled(0.)

# Get L [kg/m3]
QL = nc_file.variables['ql'][:]
L = QL*0.029*PCENTER*100./(8.314*TLAY)

# CL ICE WP
# CL ICE WP
QI = nc_file.variables['ql'][:]
REI_DEF = 24.8  #microns
RHOICE = 0.9167  #G/CM3
CICEWP = QI*(PLEV[:,:-1,:,:] - PLEV[:,1:,:,:])*100/9.8*1E3

# CL LIQ WP
RHOLIQ = 1.  # g/cm3
CLIQWP = QL*(PLEV[:,:-1,:,:] - PLEV[:,1:,:,:])*100/9.8*1E3

# REICE
REICE = np.zeros((12, LLPAR, len(lat)*len(lon))) + REI_DEF

# Read in QV
QV = nc_file.variables['qv'][:]
QV = QV *0.02896/0.0180153

QV[:,36:,:,:] = 1.E-6

#### -------------------------------------------------------------------------------####
#### Calculate liquid cloud drop radius based on CDNC change
#### -------------------------------------------------------------------------------####
# READ IN ACTIVATION
fn1 = np.load(fname1)
nact1 = fn1['CDNC']*1E6  #[m-3]
fn1.close()

fn2 = np.load(fname2)
nact2 = fn2['CDNC']*1E6  #[m-3]
fn2.close()

rdef = 10.  #[um]
RELIQ = np.zeros(shape(nact1))
mns = len(nact1)

if CONTROL:
    RELIQ[:,:,:,:] = rdef
else:
    for t in range(0, mns):
        for l in range(0, ltop):
            for la in range(0, len(lat)):
                for lo in range(0, len(lon)):
                    if PCENTER[t,l,la,lo] > 600:
                        RELIQ[t,l,la,lo] = rdef*(nact2[t,l,la,lo]/nact1[t,l,la,lo])**(1./3.)
                    else:
                        RELIQ[t,l,la,lo] = rdef
                

#### -------------------------------------------------------------------------------####
#### Make an ozone profile
#### -------------------------------------------------------------------------------####
### OZONE PROFILE TEST
ozone = np.zeros((LLPAR, len(lat), len(lon)))

# Tropopsphere ozone is 0.2 ppm
ozone[:30,:,:] = 0.2E-6

# Linearly increase o3 from trop to 10hpa
o3 = arange(0.2, 6, (6-0.2)/10.)*1E-6
ozone[30:40,:,:] = o3[:,newaxis,newaxis]

# Decrease ozone from 10hpA to top
o32 = arange(6.,1.,-(6-1)/7.)*1E-6
ozone[40:,:,:] = o32[:,newaxis,newaxis]


#### -------------------------------------------------------------------------------####
#### Reshape variables -- crop if needed, then (time, ncol)
#### -------------------------------------------------------------------------------####
# RESHAPE
t1 = 0
ltop = 30
PCENTER = PCENTER[:,:ltop,:,:]
PLEV = PLEV[:,:ltop+1,:,:]
TLAY = TLAY[:,:ltop,:,:]
TLEV = TLEV[:,:ltop+1,:,:]
TSFC = TSFC[:,:,:]
SUNCOS = SUNCOS[:,:,:,:]
clfrac = clfrac[:,:ltop,:,:]
CICEWP = CICEWP[:,:ltop,:,:]
CLIQWP = CLIQWP[:,:ltop,:,:]
RELIQ = RELIQ[:,:ltop,:,:]  
REICE = REICE[:,:ltop,:]
ozone = ozone[:ltop,:,:]
QV = QV[:,:ltop,:,:]
LLPAR = 30
time2 = 12

PCENTER = reshape(PCENTER, (time2,LLPAR, len(lat)*len(lon)))
PLEV = reshape(PLEV, (time2,LLPAR+1, len(lat)*len(lon)))
TLAY = reshape(TLAY, (time2,LLPAR,len(lat)*len(lon)))
TLEV = reshape(TLEV, (time2,LLPAR+1,len(lat)*len(lon)))
TSFC = reshape(TSFC, (time2, len(lat)*len(lon)))
ALBDRVIS = reshape(ALBDRVIS, (time2, len(lat)*len(lon)))
ALBDFVIS = reshape(ALBDFVIS, (time2, len(lat)*len(lon)))
ALBDRNIR = reshape(ALBDRNIR, (time2, len(lat)*len(lon)))
ALBDFNIR = reshape(ALBDFNIR, (time2, len(lat)*len(lon)))
CLDFR = reshape(clfrac, (time2,LLPAR, len(lat)*len(lon)))
SUNCOS = reshape(SUNCOS, (time2,HRS,len(lat)*len(lon)))
CICEWP = reshape(CICEWP, (time2, LLPAR, len(lat)*len(lon)))
CLIQWP = reshape(CLIQWP, (time2, LLPAR, len(lat)*len(lon)))
RELIQ = reshape(RELIQ, (time2,LLPAR, len(lat)*len(lon)))
O3VMR = reshape(ozone, (LLPAR, len(lat)*len(lon)))
H2OVMR = reshape(QV, (time2, LLPAR, len(lat)*len(lon)))

#### -------------------------------------------------------------------------------####
#### Write nc file
#### -------------------------------------------------------------------------------####
# Save outputs into netcdf file
nc_file = Dataset('setup_output/' + fout, 'w', format='NETCDF4')

nc_file.createDimension('month', time2)
nc_file.createDimension('wv', 30)
nc_file.createDimension('col', len(lat)*len(lon))
nc_file.createDimension('level',LLPAR)
nc_file.createDimension('levele', LLPAR+1)
nc_file.createDimension('hours', HRS)

month_nc = nc_file.createVariable('month', 'f8', ('month',))
col_nc = nc_file.createVariable('col', 'f8', ('col',))
level_nc = nc_file.createVariable('level', 'f8', ('level',))

pcenter_nc = nc_file.createVariable('pcenter', 'f8', ('month', 'level', 'col',))
plev_nc = nc_file.createVariable('plev', 'f8', ('month', 'levele', 'col',))
tsfc_nc = nc_file.createVariable('tsfc', 'f8', ('month', 'col',))
tlay_nc = nc_file.createVariable('tlay', 'f8', ('month', 'level', 'col',))
tlev_nc = nc_file.createVariable('tlev', 'f8', ('month', 'levele', 'col',))
suncos_nc = nc_file.createVariable('suncos', 'f8', ('month','hours','col'))
albdrvis_nc = nc_file.createVariable('albdrvis', 'f8', ('month', 'col',))
albdfvis_nc = nc_file.createVariable('albdfvis', 'f8', ('month', 'col',))
albdrnir_nc = nc_file.createVariable('albdrnir', 'f8', ('month', 'col',))
albdfnir_nc = nc_file.createVariable('albdfnir', 'f8', ('month', 'col',))
cldfr_nc = nc_file.createVariable('cldfr', 'f8', ('month', 'level', 'col',))
cicewp_nc = nc_file.createVariable('cicewp', 'f8', ('month', 'level', 'col',))
cliqwp_nc = nc_file.createVariable('cliqwp', 'f8', ('month', 'level', 'col',))
reice_nc = nc_file.createVariable('reice', 'f8', ('month', 'level', 'col',))
reliq_nc = nc_file.createVariable('reliq', 'f8', ('month', 'level', 'col',))
o3vmr_nc = nc_file.createVariable('o3vmr', 'f8', ('level','col'))
h2ovmr_nc = nc_file.createVariable('h2ovmr', 'f8', ('month','level','col'))

col_nc[:] = arange(0, len(lat)*len(lon),1)
level_nc[:] = arange(0, LLPAR, 1)
month_nc[:] = arange(0,time2,1)

pcenter_nc[:,:,:] = PCENTER
plev_nc[:,:,:] = PLEV
tsfc_nc[:,:] = TSFC
tlev_nc[:,:,:] = TLEV
tlay_nc[:,:,:] = TLAY
suncos_nc[:,:,:] = SUNCOS
albdrvis_nc[:,:] = ALBDRVIS
albdfvis_nc[:,:] = ALBDFVIS
albdrnir_nc[:,:] = ALBDRNIR
albdfnir_nc[:,:] = ALBDFNIR
cldfr_nc[:,:,:] = CLDFR
cicewp_nc[:,:,:] = CICEWP
cliqwp_nc[:,:,:] = CLIQWP
reice_nc[:,:,:] = REICE
reliq_nc[:,:,:] = RELIQ
o3vmr_nc[:,:] = O3VMR
h2ovmr_nc[:,:,:] = H2OVMR

nc_file.close()
