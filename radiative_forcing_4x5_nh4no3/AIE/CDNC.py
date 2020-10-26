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
directory = '../'
filename = 'ship_so2ocbcnoxscale_aw_dist.npz'

# OUTPUT FILE NAME (npz extension added in automatically)
outfile = 'aw_ship_so2ocbcnoxscale'

### UPDRAFT VELOCITY
w = 0.5 

# Model inputs
nbins = 15
nlevs = 35
nmonths = 12

# Constants
composition = ['NK','SF','SS','ECIL','ECOB','OCIL','OCOB','DUST','AW', 'NIT', 'NH4']
molwgt = [96., 58.5, 12., 12., 12., 12.,100.,18.,64., 18.,] # molecular wgt [g/mol]
kappa = [1., 1.2, 0., 0., 0.1, 0.01, 0.01, 0.67, 0.67]
#         'SF','SS','ECIL','ECOB','OCIL','OCOB','DUST','AW', 'NIT', 'NH4'
drydens = [1770.,2000.,2000.,2000.,1400.,1400.,1500., 1770., 1770.]  #kg m-3
ndry = 9
MWwater = 0.018 # mol wgt water [kg/mol]
STwater = 0.072 # surface tension of water [N m-1 or Kg s-2]
denswater = 1000. # density of water [kg m-3]
intmix = [True,True,True,False,True,True,True]  #internally mixed?
ARGbin=3  # Start parameterization at this bin

#### -------------------------------------------------------------------------------####
#### Subroutines
#### -------------------------------------------------------------------------------####

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


# This routine returns the bin limits in mass space,
# and diameters based on the assumed density
# In this script, we will calculate diameters
# from the actual mass/number ratio
# but we'll use this to get x, xk
def TOMAS_arrays(nbins=15):
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

        # Calculate mass and diameter bin edges
        if ibins == 15:
                for k in range(1,ibins+1):
                        if k < ibins-1:
                                xk[k] = xk[k-1]*4.
                        else:
                                xk[k] = xk[k-1]*32.
     
        elif ibins == 40:
                for k in range(1, ibins+1):
                        xk[k] = xk[k-1]*2
     
        x = np.sqrt(xk[:-1]*xk[1:])

        return(xk,x)


#### -------------------------------------------------------------------------------####
#### Read in TOMAS output file
#### -------------------------------------------------------------------------------####
# Read in the file 
f1 = np.load(directory + filename)
distMonthly = f1['dist'] 
lat = f1['la']
lon = f1['lo']
f1.close()

# Get mass bin definitions. Only using xk, x
# Calculate in getKappaDia now
xk, x = TOMAS_arrays(nbins=nbins)

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

ft = np.load('../met_4x5/temperature_4x5_20130101.npz') #temp met directory
temperature = ft['temperature'][:,:,:,:]
ft.close()

fs = np.load('../met_4x5/surf_press_20130101.npz') #surf press met directory
surf_press = fs['surf_press'][:,:,:]
fs.close()

# Create 3d pressure array
presse = np.zeros((nmonths,nlevs+1,len(lat),len(lon)))
for l in range(0, nlevs+1):
    presse[:,l,:,:] = ap[l] + ( bp[l] * surf_press[:,:,:])

pressure = np.zeros((nmonths,nlevs,len(lat),len(lon)))
for l in range(0, nlevs):
        pressure[:,l,:,:] = (presse[:,l] + presse[:,l+1]) / 2

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
        print(z)
        for y in range(0,len(lat)):
            for x in range(0,len(lon)):
                SMAX[m,z,y,x],CDNC2[m,z,y,x]=arg_CCN(temperature[m,z,y,x]\
                                                     , pressure[m,z,y,x]*100.
                                                     , w, dist[0,ARGbin:,z,y,x]*1E6\
                                                     , scrit[ARGbin:,z,y,x],scritl[ARGbin:,z,y,x])


#### -------------------------------------------------------------------------------####
#### SAVE OUTPUT TO FILE
#### -------------------------------------------------------------------------------####
np.savez('CDNC_output/'+outfile+'_cdnc.npz', CDNC=CDNC2, SMAX=SMAX, lat=lat, lon=lon)
