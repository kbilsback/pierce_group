# This script makes an input file for RRTMG DRE

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
# OPTICS input files -- made from 'RRTMG_coat.py" or similar
fname_optics_in = 'RRTMG_ext_output/aw_ship_so2scale_RRTMG_EXT.npz'
#fname_optics_in = 'param_elvoc_2xdms_RRTMG_COAT.npz'

# OUTFILE NAME
fout ='setup_output/aw_ship_so2scale.nc'

# Dimensions
nmonths = 12

# Cloud drop radius switch
# Unless you want to combine with AIE, leave as True
CONTROL=True

# READ IN AEROSOL
fdre = load(fname_optics_in)
SSA = fdre['SSA']
TAU = fdre['TAU']
GSCA = fdre['ASM']
fdre.close()

TAU[:,:,-1,:,:] = 0.

### RRTMG uses the reverse wavelegnth order as our scripts
SSA = SSA[:,::-1,:,:,:]
TAU = TAU[:,::-1,:,:,:]
GSCA = GSCA[:,::-1,:,:,:]

# Deal with NANS
TAU = np.ma.masked_invalid(TAU).filled(0.)
SSA = np.ma.masked_invalid(SSA).filled(1.)
GSCA = np.ma.masked_invalid(GSCA).filled(0)

#### -------------------------------------------------------------------------------####
#### Get MET variables
#### -------------------------------------------------------------------------------####
# READ IN MET FILE
nc_file = Dataset('../../met_4x5/met_albedo_4x5.nc')
lat = nc_file.variables['lat'][:]
lon = nc_file.variables['lon'][:]

LLPAR = 47

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
TLAY = TLAY[:]
TSFC = TLAY[:,0,:,:]

TLEV = np.zeros((nmonths,LLPAR+1,len(lat),len(lon)))  #FIX THIS IN THE FUTURE
TLEV[:,:-1,:,:] = TLAY[:,:,:,:]
TLEV[:,-1,:,:] = TLAY[:,-1,:,:]
TLAY = (TLEV[:,:-1,:,:]+TLEV[:,1:,:,:])/2.

# ALBEDO
#albedo_file = '/pierce-scratch/kbilsback/geos_chem/projects/china_rescoal/radiative_effects_025x03125/met_025x03125/seaice_025x03125_w1415.npz'
#nc_file.variables = np.load(albedo_file)
ALBDRVIS = nc_file.variables['albdrvis'][:]
ALBDFVIS = nc_file.variables['albdfvis'][:]
ALBDRNIR = nc_file.variables['albdrnir'][:]
ALBDFNIR = nc_file.variables['albdfnir'][:]

ALBDRVIS = ALBDRVIS[:]
ALBDFVIS = ALBDFVIS[:]
ALBDRNIR = ALBDRNIR[:]
ALBDFNIR = ALBDFNIR[:]


# CLOUD FRACTION
clfrac = nc_file.variables['clfrac'][:]
clfrac = clfrac[:]

# SOLAR DECLINATION ANGLE - LOOK INTO THIS
dhrs= 3
HRS = 24//dhrs
monSolDecAng = getMonSolDecAng()
RA = array([19.42, 21.42, 23.39, 1.32, 3.15, 5.58, 7.28, 9.41, 11.34, 13.27, 15.16, 17.4])
LONS, LATS = meshgrid(lon, lat)
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
QI = nc_file.variables['qi'][:]
REI_DEF = 24.8  #microns
RHOICE = 0.9167  #G/CM3
CICEWP = QI*(PLEV[:,:-1,:,:] - PLEV[:,1:,:,:])*100/9.8*1E3

# CL LIQ WP
RHOLIQ = 1.  # g/cm3
CLIQWP = QL*(PLEV[:,:-1,:,:] - PLEV[:,1:,:,:])*100/9.8*1E3

# REICE
REICE = zeros((12, LLPAR, len(lat),len(lon))) + REI_DEF

# Read in QV
QV = nc_file.variables['qv'][:]
QV = QV *0.02896/0.0180153

QV[:,36:,:,:] = 1.E-6

#### -------------------------------------------------------------------------------####
#### Calculate liquid cloud drop radius based on CDNC change
#### -------------------------------------------------------------------------------####
# So if you want to also do a combined AIE, read in CDNC files here
# (you can just copy the lines from rrtmg_aie_setup)

# If control is True then don't worry about it
rdef = 10.  #[um]
RELIQ = zeros((nmonths, LLPAR, len(lat), len(lon)))

if CONTROL:
    RELIQ[:,:,:,:] = rdef
else:
    print('Not coded in')
    # Copy/paste the loop from the rrtmg_aie code
    # You will need to have the CDNC files read in as well

#### -------------------------------------------------------------------------------####
#### Make an ozone profile
#### -------------------------------------------------------------------------------####
### OZONE PROFILE TEST
ozone = zeros((LLPAR, len(lat), len(lon)))

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
ltop = 30

PCENTER = PCENTER[:,:ltop,:,:]
PLEV = PLEV[:,:ltop+1,::,:]
TLAY = TLAY[:,:ltop,::,:]
TLEV = TLEV[:,:ltop+1,:,:]
TSFC = TSFC[:,:,:]
SUNCOS = SUNCOS[:,:,:,:]
clfrac = clfrac[:,:ltop,:,:]
CICEWP = CICEWP[:,:ltop,:,:]
CLIQWP = CLIQWP[:,:ltop,:,:]
RELIQ = RELIQ[:,:ltop,:,:]  
REICE = REICE[:,:ltop,:,:]
ozone = ozone[:ltop,:,:]
QV = QV[:,:ltop,:,:]

ALBDRVIS = ALBDRVIS[:,:,:]
ALBDFVIS = ALBDFVIS[:,:,:]
ALBDRNIR = ALBDRNIR[:,:,:]
ALBDFNIR = ALBDFNIR[:,:,:]

TAU = TAU[:,:,:ltop,:,:]
SSA = SSA[:,:,:ltop,:,:]
GSCA = GSCA[:,:,:ltop,:,:]

time2 = nmonths
LLPAR = 30
lat = lat[:]

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
REICE = reshape(REICE, (time2,LLPAR, len(lat)*len(lon)))
O3VMR = reshape(ozone, (LLPAR, len(lat)*len(lon)))
H2OVMR = reshape(QV, (time2, LLPAR, len(lat)*len(lon)))

TAUAER = reshape(TAU, (time2, 14, LLPAR, len(lat)*len(lon)))
#TAUAER = np.zeros(np.shape(TAUAER)) #for a simulation without any aerosol
SSAAER = reshape(SSA, (time2, 14, LLPAR, len(lat)*len(lon)))
ASMAER = reshape(GSCA, (time2, 14, LLPAR, len(lat)*len(lon)))

# Save outputs into netcdf file
nc_file = Dataset(fout, 'w', format='NETCDF4')

nc_file.createDimension('month', time2)
nc_file.createDimension('wv', 14)
nc_file.createDimension('col', len(lat)*len(lon))
nc_file.createDimension('level', LLPAR)
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

tauaer_nc = nc_file.createVariable('tauaer', 'f8', ('month','wv', 'level','col'))
ssaaer_nc = nc_file.createVariable('ssaaer', 'f8', ('month','wv', 'level','col'))
asmaer_nc = nc_file.createVariable('asmaer', 'f8', ('month','wv', 'level','col'))

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

tauaer_nc[:,:,:,:] = TAUAER
ssaaer_nc[:,:,:,:] = SSAAER
asmaer_nc[:,:,:,:] = ASMAER

nc_file.close()
