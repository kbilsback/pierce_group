# This code calculates CDNC from TOMAS
# Code calls getKappaDia and arg_CCN_func

import numpy as np
from netCDF4 import Dataset
from getKappaDia import getKappaDia
from arg_CCN_func import arg_CCN
import sys

#### -------------------------------------------------------------------------------####
#### Constants and inputs
#### -------------------------------------------------------------------------------####
# path to file
directory = '../test_simulations/'
#filename = 'param_ELVOC_2xDMS_dist_monthly.npz'
filename = 'base_noMSA_dist_monthly.npz'

# OUTPUT FILE NAME (npz extension added in automatically)
#outfile = 'param_elvoc_2xdms'
outfile = 'base'

### UPDRAFT VELOCITY
w = 0.5 

# Model inputs
nbins = 15
nlevs = 35
nmonths = 12

# Constants
composition = ['NK','SF','SS','ECIL','ECOB','OCIL','OCOB','DUST','AW']
molwgt = [96., 58.5, 12., 12., 12., 12.,100.,18.] # molecular wgt [g/mol]
kappa = [1.,1.2,0.,0.,0.1,0.01,0.01]
drydens = [1770.,2000.,2000.,2000.,1400.,1400.,1500.]  #kg m-3
ndry = 7
MWwater = 0.018 # mol wgt water [kg/mol]
STwater = 0.072 # surface tension of water [N m-1 or Kg s-2]
denswater = 1000. # density of water [kg m-3]
intmix = [True,True,True,False,True,True,True]  #internally mixed?
ARGbin=3  # Start parameterization at this bin

#### -------------------------------------------------------------------------------####
#### Subroutines
#### -------------------------------------------------------------------------------####
# Pressure middles (average/standard)
def get_presm():
        etam=np.array([0.000028,0.000127,0.000400,0.001113,0.002825,0.006609,\
                                          0.014389,0.023832,0.033923,0.047794,0.066773,0.084584,\
                                          0.099510,0.117070,0.137729,0.162033,0.190668,0.224471,\
                                          0.264374,0.310730,0.354300,0.391937,0.429595,0.467274,\
                                          0.504968,0.542673,0.580389,0.618113,0.655841,0.687287,\
                                          0.712447,0.737606,0.762768,0.787933,0.810582,0.828198,\
                                          0.843298,0.858399,0.873500,0.888601,0.903703,0.918805,\
                                          0.933908,0.949010,0.964113,0.979217,0.994283]) # box middle

        ptop=0.01 #mb
        psurf=1000 #mb
        nlevs = 47
        presm = []
        for ll in range(0,nlevs):
                l=nlevs-ll-1
                presm.append(ptop+etam[l]*(psurf-ptop))

        presm=np.array(presm)
        return presm

# Pressure edges (average/standard)
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

# This routine returns the bin limits in mass space,
# and diameters based on the assumed density
# In this script, we will calculate diameters
# from the actual mass/number ratio
# but we'll use this to get x, xk
def TOMAS_arrays(nbins=15, dens=1400.):
        # Wrriten by R. Stevens
        # Returns:
        #   xk: bin edges in mass space [kg]
        #   x: average mass in each bin [kg]
        #   Dpk: bin edges in diameter [um]
        #   Dp: average diameter of each bin [um]

        ibins = nbins
        #    dens = 2160.   sea salt
        #    dens = 1700.
        #    dens = 2200. # assuming (NH4)2SO4 [kg/m3]
        #    dens = 1400.  #test

        xk = np.zeros(ibins+1)  # particle mass cutoffs [kg]
        Dpk = np.zeros(ibins+1) # particle diameter cutoffs [um]

        if ibins ==15:
                xk[0]=1.6E-23  # avg mass per particle in lowest bin edge
        elif ibins == 40:
                xk[0] = 1.0E-21 * 2.0**(-10)

        vol = xk[0]/dens
        Dpk[0]=(6.*vol/np.pi)**(1./3)*1e6

        # Calculate mass and diameter bin edges
        if ibins == 15:
                for k in range(1,ibins+1):
                        if k < ibins-1:
                                xk[k] = xk[k-1]*4.
                        else:
                                xk[k] = xk[k-1]*32.
                                vol = xk[k]/dens
                                Dpk[k]=(6.*vol/np.pi)**(1./3)*1e6

        elif ibins == 40:
                for k in range(1, ibins+1):
                        xk[k] = xk[k-1]*2
                        vol = xk[k]/dens
                        Dpk[k]=(6.*vol/np.pi)**(1./3)*1e6

        x = np.sqrt(xk[:-1]*xk[1:])
        vol = x/dens
        Dp = (6.*vol/np.pi)**(1./3)*1e6

        return(xk,x,Dpk,Dp)

#### -------------------------------------------------------------------------------####
#### Read in TOMAS output file
#### -------------------------------------------------------------------------------####
# Read in the file 
f1 = np.load(directory + filename)
distMonthly = f1['dist'] 
lat = f1['lat']
lon = f1['lon']
f1.close()

# Get mass bin definitions. Only using xk, x
xk, x, dpkF, dpF = TOMAS_arrays(nbins=nbins)

#### If you haven't already:
#### Change units here to num cm-3 and ug m-3
### Volume should be for STP (convert to ambient later)
### Or comment out line later that converts to ambient
# Number species
distMonthly[:,0,:,:,:,:] = distMonthly[:,0,:,:,:,:]/22.4e15*273./293.   #num cm-3 at STP

# Mass species
for c in range(0,len(molwgt)):
        distMonthly[:,c+1,:,:,:,:] = distMonthly[:,c+1,:,:,:,:]*molwgt[c]/22.4e15*273/293.*1e6*1E9  #ug m-3

#### -------------------------------------------------------------------------------####
#### Read in temperature and pressure files
#### -------------------------------------------------------------------------------####

ft = np.load('/home/jkodros/gloo/gregg/met/temperature2010_4x5.npz') #temp met directory
temperature = ft['temperature'][:,:,:,:]
ft.close()

fs = np.load('/home/jkodros/gloo/gregg/met/surf_press2010_4x5.npz') #surf press met directory
surf_press = fs['surf_press'][:,:,:]
fs.close()

presm = get_presm()
# Create 3d pressure array
pressure = np.zeros((nmonths,nlevs,len(lat),len(lon)))
for l in range(0, nlevs):
    pressure[:,l,:,:] = surf_press[:,:,:]*presm[l]/presm[0]

#### -------------------------------------------------------------------------------####
#### Crop files to dimensions of interest
#### -------------------------------------------------------------------------------####
# You may wish to do the calculation only on a subset of the grid
temperature = temperature[:,:nlevs,:,:]
pressure = pressure[:,:nlevs,:,:]
distMonthly = distMonthly[:,:,:,:nlevs,:,:]

#### -------------------------------------------------------------------------------####
#### MAIN CODE LOOP
#### -------------------------------------------------------------------------------####

# Calc A factor
Afactor = np.zeros(np.shape(temperature))
Afactor=4.*STwater*MWwater/(8.314*temperature*denswater)

# Create CDNC and SMAX output arrays
CDNC2 = np.zeros((nmonths, nlevs, len(lat), len(lon)))
SMAX = np.zeros((nmonths, nlevs, len(lat), len(lon)))

# Loop over months
for m in range(0, nmonths):
    dist = distMonthly[m,:,:,:,:,:]

    # Convert from STP to ambient concentrations
    dist = dist[:,:,:,:,:]*273/temperature[np.newaxis,m,:,:,:]*pressure[np.newaxis,m,:,:,:]/1E3

    # Change mass components to kg cm3
    dist[1:,:,:,:,:] = dist[1:,:,:,:,:]*1E-9*1E-6

    # Call the getKappDia code
    # This calculates average kappa and bin diamaters
    ddia,ddialo,ddiahi,tkappa,dmass = getKappaDia(dist,kappa,drydens,xk,nbins,ndry)

    #### Calculate Scrit
    scrit=np.exp(np.sqrt(4.*Afactor[m,np.newaxis,:,:,:]**3/(27.*ddia**3*tkappa)))-1.
    # scrit at limits
    scritl=np.empty((nbins+1,nlevs,len(lat),len(lon)))
    scritl[0,:,:,:]=(np.exp(np.sqrt(4.*Afactor[m,:,:,:]**3/(27.*ddialo[0,:,:,:]**3*tkappa[0,:,:,:])))-1.)
    scritl[1:-1,:,:,:]=(np.exp(np.sqrt(4.*Afactor[m,:,:,:]**3/(27.*((ddiahi[:-1,:,:,:]+ddialo[1:,:,:,:])/2.)**3*((tkappa[:-1,:,:,:]+tkappa[1:,:,:,:])/2.))))-1.)
    scritl[-1,:,:,:]=(np.exp(np.sqrt(4.*Afactor[m,:,:,:]**3/(27.*ddiahi[-1,:,:,:]**3*tkappa[-1,:,:,:])))-1.)

    # Loop over lev,lat,lon and calculate CDNC and SMAX
    for z in range(0,nlevs):
        # for time keeping purposes, print level
        print z
        for y in range(0,len(lat)):
            for x in range(0,len(lon)):
                SMAX[m,z,y,x],CDNC2[m,z,y,x]=arg_CCN(temperature[m,z,y,x]\
                                                     , pressure[m,z,y,x]*100.
                                                     , w, dist[0,ARGbin:,z,y,x]*1E6\
                                                     , scrit[ARGbin:,z,y,x],scritl[ARGbin:,z,y,x])


#### -------------------------------------------------------------------------------####
#### SAVE OUTPUT TO FILE
#### -------------------------------------------------------------------------------####
np.savez('CDNC_output/'outfile+'_cdnc.npz', CDNC=CDNC2, SMAX=SMAX, lat=lat, lon=lon)
