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

def get_prese():
    etae=np.array([0.000000,0.000055,0.000199,0.000601,0.001625,0.004026,\
                   0.009191,0.019586,0.028077,0.039768,0.055820,0.077726,\
                   0.091442,0.107578,0.126563,0.148896,0.175170,0.206167,\
                   0.242774,0.285974,0.335486,0.373114,0.410759,0.448431,\
                   0.486118,0.523819,0.561527,0.599251,0.636974,0.674708,\
                   0.699867,0.725026,0.750186,0.775350,0.800515,0.820648,\
                   0.835748,0.850848,0.865949,0.881051,0.896152,0.911253,\
                   0.926356,0.941459,0.956562,0.971665,0.986769,1.001796])# box edge sigma coordinates

    ptop=0.01 #mb
    psurf=1000 #mb
    nlevs=47
    prese=[]
    for ll in range(0,nlevs+1):
        l=(nlevs+1)-ll-1
        prese.append(ptop+etae[l]*(psurf-ptop))
    prese=np.array(prese)
    return prese

#### -------------------------------------------------------------------------------####
#### File names and options
#### -------------------------------------------------------------------------------####
# CDNC FILES -- made from CDNC.py
fname1 = 'param_elvoc_2xdms_cdnc.npz'
fname2 = 'base_cdnc.npz'

# OUTFILE NAME
fout ='rrtmg_inputs_CTRL.nc'
CONTROL=True

# Only do rad-transfer on first 30 levels
ltop = 30
#### -------------------------------------------------------------------------------####
#### Get MET variables
#### -------------------------------------------------------------------------------####
# READ IN MET FILE
nc_file = Dataset('../met_4x5/met_4x5.nc')
lat = nc_file.variables['lat'][:]
lon = nc_file.variables['lon'][:]

LLPAR = 47

# Pressure Endge
surf_press = nc_file.variables['surf_press'][:]  #[time,lat,lon]
prese = get_prese()
PEDGE = np.zeros((12,LLPAR+1,len(lat),len(lon)))
PLEV = np.zeros((12,LLPAR+1,len(lat),len(lon)))
for i in range(0, len(prese)):
    PEDGE[:,i,:,:] = surf_press[:,:,:] * prese[i]/prese[0]
PLEV[:,:,:,:] = PEDGE[:,:,:,:]
PCENTER = (PLEV[:,:-1,:,:]+PLEV[:,1:,:,:])/2.

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
HRS = 24/dhrs
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
QL = nc_file.variables['QL'][:]
L = QL*0.029*PCENTER*100./(8.314*TLAY)

# CL ICE WP
# CL ICE WP
QI = nc_file.variables['QI'][:]
REI_DEF = 24.8  #microns
RHOICE = 0.9167  #G/CM3
CICEWP = QI*(PLEV[:,:-1,:,:] - PLEV[:,1:,:,:])*100/9.8*1E3

# CL LIQ WP
RHOLIQ = 1.  # g/cm3
CLIQWP = QL*(PLEV[:,:-1,:,:] - PLEV[:,1:,:,:])*100/9.8*1E3

# REICE
REICE = np.zeros((12, LLPAR, len(lat)*len(lon))) + REI_DEF

# Read in QV
QV = nc_file.variables['QV'][:]
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
