# This script calculates AOD, SSA, ASYM for TOMAS output
# This is intended to be used with RRTMG to calc DRE offline
# Core-shell assumption for BC

# Written by JKodros 9/4/14
# Updated JKodros 10/6/14 - BHCOAT now returns g
# Updated JKodros 5/19/16 - Module use on ozone
# Updated JKodros 9/26/17 - Re-writing to make easier for others to use

import os
import numpy as np
import bhcoat_callg
import sys

#### -------------------------------------------------------------------------------####
#### Constants and inputs
#### -------------------------------------------------------------------------------####
# path to file
directory = '/pierce-scratch/kbilsback/geos_chem/projects/china_rescoal/radiative_effects_025x03125/'
fname = 'china_boff_hi_dist_w1415.npz'

# OUTPUT FILE NAME (npz extension added in automatically)
outf = 'china_boff_hi_w1415'

### OPTIONS
BROWN=False

# Model inputs
nbins = 15
nlevs = 47
nmonths = 12

# Constants
composition = ['NK','SF','SS','ECIL','ECOB','OCIL','OCOB','DUST','AW']
molwgt = [96., 58.5, 12., 12., 12., 12.,100.,18.] # molecular wgt [g/mol]
kappa = [1.,1.2,0.,0.,0.1,0.01,0.01]
ndry = 7
grav=9.814 # Gravity [m s-2]

# Aerosol Density
rho_so4=1.78e3  #I.N. Tang (1996) [kg/m3] (Ammonium bisulfate = 1780kg/m3)
rho_h2o=1.e3
rho_oc=1.4e3  # Dick et al. (2000) - consistent with water uptake reference
rho_bc=1.8e3  # Bond and Bergstrom (2006)
rho_dust=2.65e3 #  Tegen and Fung (1994)
rho_nacl=2.165e3 # I.N. Tang (1996)

# SW Wavelengths used in RRTMG
wl = [229.8, 329.1,399.8,527.1,693.5,944.3,1270.3,1455.2,1778.4,2044.2,2320.2,2777.0,3418.8, 3444.7]

# Read in refractive indicies
frf = np.load('ref_ind.npz')
n_file = frf['n']
k_file = frf['k']
frf.close()

# Assign RFs to variables
n_so4 = np.array(n_file[0][:14])
n_nacl= np.array(n_file[1][:14])
n_bc =  np.array(n_file[2][:14])
n_oc =  np.array(n_file[3][:14])
n_dust= np.array(n_file[4][:14])
n_h2o = np.array(n_file[5][:14])

k_so4 = abs(np.array(k_file[0][:14]))
k_nacl= abs(np.array(k_file[1][:14]))
k_bc =  abs(np.array(k_file[2][:14]))
k_oc =  abs(np.array(k_file[3][:14]))
k_dust= abs(np.array(k_file[4][:14]))
k_h2o = abs(np.array(k_file[5][:14]))

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
#### MAIN CODE
#### -------------------------------------------------------------------------------####
# Read in TOMAS output file
# KEEP IN ppb (unless you edit the code below yourself)
f1 = np.load(directory + fname)
dist = f1['dist']
lat = f1['la']
lon = f1['lo']
f1.close()

# DEFINE OUTPUT ARRAYS
TAU = np.zeros((nmonths,len(wl),nlevs,len(lat),len(lon)))
SSA = np.zeros((nmonths,len(wl),nlevs,len(lat),len(lon)))
ASM = np.zeros((nmonths,len(wl),nlevs,len(lat),len(lon)))

print('Starting DRE Calculation')
# Main time loop
for t in range(0, nmonths):
    print('Month ' + str(t))

    ### Get total aerosol number in num (kg air)-1 ###
    # Mask very low values
    aeronum = dist[t,0,:,:,:,:]/22.4e15*1e6/1.28   # num/kg air
    aeronum = np.ma.masked_less(aeronum, 1E-10)
    aeronum = np.ma.filled(aeronum, 0.)

    ### Get mass species in kg/kg air ###
    # Put into one list "mass"
    for c in range(0,len(molwgt)):
        dist[t,c+1,:,:,:,:] = dist[t,c+1,:,:,:,:]*molwgt[c]/22.4e15*1e6/1.28   #kg/kgair

    mso4=np.ma.masked_array(dist[t,1,:,:,:,:])
    mnh4=mso4*0.1875
    mnacl=np.ma.masked_array(dist[t,2,:,:,:,:])
    mecil=np.ma.masked_array(dist[t,3,:,:,:,:])
    mecob=np.ma.masked_array(dist[t,4,:,:,:,:])
    moc=np.ma.masked_array(dist[t,5,:,:,:,:]+dist[t,6,:,:,:,:])
    mdust=np.ma.masked_array(dist[t,7,:,:,:,:])
    mh2o=np.ma.masked_array(dist[t,8,:,:,:,:])
    mtot=(mso4+mnh4+mnacl+moc+mecil+mecob+mdust) # Total dry mass
    mp=mtot/aeronum # Mass per particle
    #mass = [mso4, mnh4, mnacl, mecil, mecob, moc, mdust, mh2o, mtot,  mp]

    ### Get volumes ###
    # VOLUMES m3
    vol_so4  = (mso4+mnh4)/rho_so4 #m3
    vol_h2o = mh2o/rho_h2o
    vol_oc  = moc/rho_oc
    vol_dust= mdust/rho_dust
    vol_nacl= mnacl/rho_nacl
    vol_ecob=mecob/rho_bc
    vol_ecil=mecil/rho_bc
    #vol = [vol_so4, vol_h2o, vol_oc, vol_dust, vol_nacl, vol_ecob, vol_ecil]

    ### LOOP OVER WAVELENGTHS
    for w in range(0, len(wl)):
        print('wl ' + str(wl[w]))

        rfcore = complex(n_bc[w], k_bc[w])

        #### Get refractive indicies for this wavelength ###
        # If brC, use Saleh et al. (2014) for OA
        if BROWN:
            mass_ecil = np.array(mass[3])
            mass_ecob = np.array(mass[4])
            mass_oc = np.array(mass[5])
        
            BCOA_ratio = np.zeros(np.shape(mass_ecil))
        
            BCOA_ratio[:,:,:,:] = (mass_ecil[:,:,:,:]+mass_ecob[:,:,:,:])/(mass_oc[:,:,:,:])
            k = rad.calc_kOA(BCOA_ratio, wl[w])
            k2 = np.ma.masked_invalid(k)
            k2 = np.ma.masked_greater(k2, 0.05)
            k2 = np.ma.filled(k2, 0.05)
            k3 = np.ma.masked_less(k2, 6.E-3)
            k_oc_new = np.ma.filled(k3, 6.E-3)
        
        ### Get refractive index of the scattering species ###
        # i.e. the mantle/shell
        
        # Volume of external mixture
        vol_ext = vol_ecob + vol_ecil

        #Volume of internal
        vol_int = vol_so4+vol_oc+vol_dust+vol_nacl

        # Real index of refraction
        # NOTE: we are including aerosol water in this calc
        # NOTE: for brown carbon replace k_oc with k_oc_new (calculated above)
        refre_mant=(vol_so4*n_so4[w]+vol_h2o*n_h2o[w]+vol_nacl*n_nacl[w]+vol_oc*n_oc[w]\
                    +vol_dust*n_dust[w])/(vol_int+vol_h2o)
        # Imaginary Refractive Index
        refim_mant=(vol_so4*k_so4[w]+vol_h2o*k_h2o[w]+vol_nacl*k_nacl[w]+vol_oc*k_oc[w]\
                    +vol_dust*k_dust[w])/(vol_int+vol_h2o)

        # Complex number
        rfmant = refre_mant + refim_mant*1j

        ### Calculate radius of the core and shell ###
        # BC (core) volume per particle 
        vol_corep = (vol_ecob + vol_ecil)/aeronum
        r_core = (3./4.*vol_corep/np.pi)**(1./3.)

        # Shell volume per particle (con agua)
        vol_totp = (vol_so4+vol_oc+vol_dust+vol_nacl+vol_ecob+vol_ecil+vol_h2o)/aeronum
        r_tot = (3./4.*vol_totp/np.pi)**(1./3.)
        r_mant = r_tot - r_core

        ### Calculate size parameter ###
        # size parameter is unitless
        sizep_core = 2*np.pi*r_core/(wl[w]*1e-9)
        sizep_mant = 2*np.pi*r_mant/(wl[w]*1e-9)
        sizep_tot = 2*np.pi*r_tot/(wl[w]*1e-9)

        # Replace masked values with 0
        sizep_core = np.ma.masked_invalid(sizep_core)
        sizep_core = np.ma.filled(sizep_core, 0.)
        sizep_mant = np.ma.masked_invalid(sizep_mant)
        sizep_mant = np.ma.filled(sizep_mant, 0.)
        sizep_tot = np.ma.masked_invalid(sizep_tot)
        sizep_tot = np.ma.filled(sizep_tot, 0.)

        ### CALL TO BHCOAT FORTRAN CODE ###
        Qbhcoat = bhcoat_callg.bhcoatloop(sizep_core, sizep_tot, rfcore, rfmant)
        qext = np.array(Qbhcoat[0])         # extinction ratio
        qsca = np.array(Qbhcoat[1])         # scattering ratio
        gsca_int = np.array(Qbhcoat[3])     # Assymetry param
        qabs = np.zeros(np.shape(qext))      
        qabs[:,:,:,:] = qext[:,:,:,:] - qsca[:,:,:,:]  # absorption ratio

        ### Calculate abs/scat coefficient ###
        b_abs = (r_tot)**2*np.pi*aeronum*qabs
        b_sca = (r_tot)**2*np.pi*aeronum*qsca

        ### Calc abs/sca AOD
        # Get average pressure levels
        presm = get_presm()
        pt = presm[np.newaxis,:,np.newaxis,np.newaxis]
        prese = get_prese()
        pet = prese[np.newaxis,:,np.newaxis,np.newaxis]

        tau_abs = b_abs*pt/grav*100.*np.log(pet[:,:-1,:,:]/pet[:,1:,:,:])
        tau_sca = b_sca*pt/grav*100.*np.log(pet[:,:-1,:,:]/pet[:,1:,:,:])

        ### Pass AOD, SSA, ASYM to output arrays ###
        TAU[t,w,:,:,:] = tau_sca[:,:,:,:].sum(0) + tau_abs[:,:,:,:].sum(0)
        SSA[t,w,:,:,:] = tau_sca[:,:,:,:].sum(0)/(tau_sca[:,:,:,:].sum(0) + tau_abs[:,:,:,:].sum(0))
        ASM[t,w,:,:,:] = (gsca_int[:,:,:,:]*b_sca[:,:,:,:]).sum(0)/b_sca.sum(0)

# Save output
np.savez('/pierce-scratch/kbilsback/geos_chem/projects/china_rescoal/radiative_effects_025x03125/DRE/core-shell/RRTMG_coat_output/'+outf+'_RRTMG_COAT.npz',lat=lat,lon=lon, TAU=TAU, SSA=SSA, ASM=ASM)
