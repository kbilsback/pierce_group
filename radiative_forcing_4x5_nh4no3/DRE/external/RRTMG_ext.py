# This script calculates AOD, SSA, ASYM for TOMAS output
# This is intended to be used with RRTMG to calc DRE offline
# External assumption for BC

# Uses bhmie

# Written by JKodros 9/4/14 
# UPDATED 10/7/14 to have browning
# Updated JKodros 4/20/18 - Re-writing to make easier for others to use

import numpy as np
import sys
import bhloop2 as bh
#import bhloop2

#### -------------------------------------------------------------------------------####
#### Constants and inputs
#### -------------------------------------------------------------------------------####
# path to file
directory = '/pierce-scratch/kbilsback/geos_chem/projects/shipping_paper/radiative_forcing_4x5_nit/'
fname = 'ship_noscale_aw_dist.npz'

# OUTPUT FILE NAME (npz extension added in automatically)
outf = 'aw_ship_noscale'

### OPTIONS
BROWN=False
ENHANCE=False

# Model inputs
nbins = 15
nlevs = 47
nmonths = 12

# Constants
composition = ['NK','SF','SS','ECIL','ECOB','OCIL','OCOB','DUST','AW', 'NIT' , 'NH4']
molwgt = [96., 58.5, 12., 12., 12., 12., 100., 18., 62., 18.] # molecular wgt [g/mol]
kappa = [1., 1.2, 0., 0., 0.1, 0.01, 0.01]
ndry = 7
grav=9.814 # Gravity [m s-2]

# Aerosol Density
rho_so4=1.78e3  #I.N. Tang (1996) [kg/m3] (Ammonium bisulfate = 1780kg/m3)
rho_h2o=1.e3
rho_oc=1.4e3  # Dick et al. (2000) - consistent with water uptake reference
rho_bc=1.8e3  # Bond and Bergstrom (2006)
rho_dust=2.65e3 #  Tegen and Fung (1994)
rho_nacl=2.165e3 # I.N. Tang (1996)
rho_nit=rho_so4
rho_nh4=rho_so4

# SW Wavelengths used in RRTMG
wl = [229.8, 329.1,399.8,527.1,693.5,944.3,1270.3,1455.2,1778.4,2044.2,2320.2,2777.0,3418.8, 3444.7]

# Read in refractive indicies
frf = np.load('ref_ind.npz')
n_file = frf['n']
k_file = frf['k']

# Assign RFs to variables
n_so4 = np.array(n_file[0][:14])
n_nacl= np.array(n_file[1][:14])
n_bc =  np.array(n_file[2][:14])
n_oc =  np.array(n_file[3][:14])
n_dust= np.array(n_file[4][:14])
n_h2o = np.array(n_file[5][:14])
n_nit = n_so4
n_nh4 = n_so4

k_so4 = abs(np.array(k_file[0][:14]))
k_nacl= abs(np.array(k_file[1][:14]))
k_bc =  abs(np.array(k_file[2][:14]))
k_oc =  abs(np.array(k_file[3][:14]))
k_dust= abs(np.array(k_file[4][:14]))
k_h2o = abs(np.array(k_file[5][:14]))
k_nit = k_so4
k_nh4 =k_so4

#### -------------------------------------------------------------------------------####
#### Subroutines
#### -------------------------------------------------------------------------------####
fs = np.load('../../met_4x5/surf_press_20130101.npz') #surf press met directory
surf_press = fs['surf_press'][:,:,:]
fs.close()

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
   for c in range(np.shape(dist)[1]-1): #range(0,len(molwgt)):
      dist[t,c+1,:,:,:,:] = dist[t,c+1,:,:,:,:]*molwgt[c]/22.4e15*1e6/1.28   #kg/kgair


   mso4=np.ma.masked_array(dist[t,1,:,:,:,:])
   #mnh4=mso4*0.1875
   mnacl=np.ma.masked_array(dist[t,2,:,:,:,:])
   mecil=np.ma.masked_array(dist[t,3,:,:,:,:])
   mecob=np.ma.masked_array(dist[t,4,:,:,:,:])
   moc=np.ma.masked_array(dist[t,5,:,:,:,:]+dist[t,6,:,:,:,:])
   mdust=np.ma.masked_array(dist[t,7,:,:,:,:])
   mh2o=np.ma.masked_array(dist[t,8,:,:,:,:])
   mnit=np.ma.masked_array(dist[t,9,:,:,:,:])
   mnh4=np.ma.masked_array(dist[t,10,:,:,:,:])
   mtot=(mso4 + mnacl + moc + mecil + mecob + mdust) # Total dry mass
   mp=mtot/aeronum # Mass per particle

   ### Get volumes ###
   # VOLUMES m3
   vol_so4  = mso4/rho_so4 #m3
   vol_h2o = mh2o/rho_h2o
   vol_nit = mnit/rho_nit
   vol_nh4 = mnh4/rho_nh4
   vol_oc  = moc/rho_oc
   vol_dust = mdust/rho_dust
   vol_nacl = mnacl/rho_nacl
   vol_ecob = mecob/rho_bc
   vol_ecil = mecil/rho_bc

   # Calc number of BC and internal populations
   num_ext = mecil[:,:,:,:]/mp[:,:,:,:] + mecob[:, :,:,:]/mp[:,:,:,:]
   num_int = aeronum[:,:,:,:]-num_ext[:,:,:,:]
   num_int = np.ma.masked_less(num_int,1E-10).filled(0.)
   num_ext = np.ma.masked_less(num_ext,1E-10).filled(0.)
            
   ### LOOP OVER WAVELENGTHS
   for w in range(0, len(wl)):
      print('wl ' + str(wl[w]))

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
         k3 = np.ma.masked_less(k2, 1.E-3)
         k_oc_new = np.ma.filled(k3, 6.E-3)

      ### Get refractive index of the scattering species ###
      # i.e. the mantle/shell

      # Volume of external mixture
      vol_ext = vol_ecob + vol_ecil

      #Volume of internal
      vol_int = vol_so4 + vol_oc + vol_dust + vol_nacl

      # Real index of refraction
      # NOTE: we are including aerosol water in this calc
      # NOTE: for brown carbon replace k_oc with k_oc_new (calculated above)
      refre = (vol_so4*n_so4[w]+vol_h2o*n_h2o[w]+vol_nacl*n_nacl[w]+vol_oc*n_oc[w]\
                  +vol_dust*n_dust[w]+vol_nit*n_nit[w]+vol_nh4*n_nh4[w])/(vol_int+vol_h2o+vol_nit+vol_nh4)
      # Imaginary Refractive Index
      refim = (vol_so4*k_so4[w]+vol_h2o*k_h2o[w]+vol_nacl*k_nacl[w]+vol_oc*k_oc[w]\
                  +vol_dust*k_dust[w]+vol_nit*k_nit[w]+vol_nh4*k_nh4[w])/(vol_int+vol_h2o+vol_nit+vol_nh4)

      # Complex number
      rfint = refre + refim*1j

      # Calculate radius of internal and external populations
      vol_extp = (vol_ecil+vol_ecob)/num_ext
      r_ext = (3./4.*vol_extp/np.pi)**(1./3.)

      vol_intp = (vol_so4+vol_h2o+vol_oc+vol_dust+vol_nacl+vol_nh4+vol_nit)/num_int
      r_int = (3./4.*vol_intp/np.pi)**(1./3.)

      ### Calculate size parameter ###
      # size parameter is unitless
      sizep_int = 2*np.pi*r_int/(wl[w]*1e-9)
      sizep_ext = 2*np.pi*r_ext/(wl[w]*1e-9)

      # Replace masked values with 0
      sizep_int = np.ma.masked_invalid(sizep_int)
      sizep_int = np.ma.filled(sizep_int, 0.)
      sizep_ext = np.ma.masked_invalid(sizep_ext)
      sizep_ext = np.ma.filled(sizep_ext, 0.)

      # Comples rfracive index (just that of BC)
      rfre_ext = np.zeros(np.shape(refre)) + n_bc[w]
      rfim_ext = np.zeros(np.shape(refre)) + k_bc[w]
      rfext = np.vectorize(complex)(rfre_ext, rfim_ext)

      ### CALL TO BHCOAT FORTRAN CODE ###
      # Internal population
      Qint = bh.bhloop(sizep_int, rfint)
      qext_int = np.array(Qint[0])
      qsca_int = np.array(Qint[1])
      gsca_int = np.array(Qint[2])
      qabs_int = np.zeros(np.shape(qext_int))
      qabs_int[:,:,:,:] = qext_int[:,:,:,:] - qsca_int[:,:,:,:]

      # External Population
      # (note its confusting 'ext' is either external or extiction
      Qext = bh.bhloop(sizep_ext, rfext)
      qext_ext = np.array(Qext[0])
      qsca_ext = np.array(Qext[1])
      gsca_ext = np.array(Qext[2])
      qabs_ext = np.zeros(np.shape(qext_ext))
      qabs_ext[:,:,:,:] = qext_ext[:,:,:,:] - qsca_ext[:,:,:,:]
      if ENHANCE:
         qabs_ext = qabs_ext * 1.5

      #### Get average pressure levels
      ### (this is for aod)
      #print('Here!')
      # Create 3d pressure array
      prese = np.zeros((nlevs+1, len(lat), len(lon)))
      for l in range(0, nlevs+1):
         prese[l,:,:] = ap[l] + ( bp[l] * surf_press[t,:,:])

      presm = np.zeros((nlevs, len(lat), len(lon)))
      for l in range(0, nlevs):
         presm[l,:,:] = (prese[l] + prese[l+1]) / 2
      pt = presm[np.newaxis,:]
      pet = prese[np.newaxis,:]

      ### INTERNAL POPULATION
      b_abs_int = (r_int)**2*np.pi*num_int*qabs_int
      b_sca_int = (r_int)**2*np.pi*num_int*qsca_int

      tau_sca_int = b_sca_int*pt/grav*100.*np.log(pet[:,:-1,:,:]/pet[:,1:,:,:])
      tau_abs_int = b_abs_int*pt/grav*100.*np.log(pet[:,:-1,:,:]/pet[:,1:,:,:])

      tau_int = tau_sca_int[:,:,:,:].sum(0) + tau_abs_int[:,:,:,:].sum(0)

      ssa_int = tau_sca_int[:,:,:,:].sum(0)/(tau_int[:,:,:])

      gsca_int_all = (gsca_int[:,:,:,:]*b_sca_int[:,:,:,:]).sum(0)/b_sca_int.sum(0)

      ### EXTERNAL POPULATION
      b_abs_ext = (r_ext)**2*np.pi*num_ext*qabs_ext
      b_sca_ext = (r_ext)**2*np.pi*num_ext*qsca_ext

      tau_sca_ext = b_sca_ext*pt/grav*100.*np.log(pet[:,:-1,:,:]/pet[:,1:,:,:])
      tau_abs_ext = b_abs_ext*pt/grav*100.*np.log(pet[:,:-1,:,:]/pet[:,1:,:,:])

      tau_ext = tau_sca_ext[:,:,:,:].sum(0) + tau_abs_ext[:,:,:,:].sum(0)

      ssa_ext = tau_sca_ext[:,:,:,:].sum(0)/(tau_ext[:,:,:])

      gsca_ext_all = (gsca_ext[:,:,:,:]*b_sca_ext[:,:,:,:]).sum(0)/b_sca_ext.sum(0)

      # PUTTING IT TOGETHER
      # see klingmuller et al 2014
      TAU[t,w,:,:,:] = tau_int[:,:,:] + tau_ext[:,:,:]

      SSA[t,w,:,:,:] = (tau_int[:,:,:]*ssa_int[:,:,:] + tau_ext[:,:,:]*ssa_ext[:,:,:])\
          /(tau_int[:,:,:] + tau_ext[:,:,:])

      ASM[t,w,:,:,:] = (tau_int[:,:,:]*ssa_int[:,:,:]*gsca_int_all[:,:,:] + tau_ext[:,:,:]*ssa_ext[:,:,:]*gsca_ext_all[:,:,:])\
          /(tau_int[:,:,:]*ssa_int[:,:,:] + tau_ext[:,:,:]*ssa_ext[:,:,:])

# Save output
np.savez('RRTMG_ext_output/'+outf+'_RRTMG_EXT.npz',lat=lat,lon=lon, TAU=TAU, SSA=SSA, ASM=ASM)



