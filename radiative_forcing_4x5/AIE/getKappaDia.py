# getKappaDia: This function takes the number and mass distributions and
# calculates the effective hygroscopisity parameters (kappa) and the dry
# diameter for particles at the low and high limit of each size bin

import numpy as np
#from numba import autojit

#@autojit
def getKappaDia(dist,kappa,drydens,xk,nbins,ndry):

   Nko=dist[0,:,:,:,:]   #number
   Mko=dist[1:,:,:,:,:]  #mass
   
   pvol=[] # partial volume per component
   for c in range(0,ndry):
      pvol.append(Mko[c,:,:,:,:]/Nko/drydens[c])

   pvol=np.array(pvol)
#   print 'pvol.shape',pvol.shape  
   dvol=pvol.sum(0) # total dry volume
   ddia=(6.0*dvol/np.pi)**(1.0/3.0) # dry diameter
   
   dmass = Mko[:ndry,:,:,:,:].sum(0)/Nko # dry mass per particle
   
   ddialo=np.empty(ddia.shape)
   ddiahi=np.empty(ddia.shape)
   for k in range(0,nbins): 
      ddialo[k,:,:,:] = ddia[k,:,:,:]*(xk[k]/dmass[k,:,:,:])**(1.0/3.0)
      ddiahi[k,:,:,:] = ddia[k,:,:,:]*(xk[k+1]/dmass[k,:,:,:])**(1.0/3.0)
 #  print 'ddialo.shape',ddialo.shape  
   
   tkappa=np.zeros(Nko.shape)
   for c in range(0,ndry):
      tkappa=tkappa+pvol[c,:,:,:,:]/dvol*kappa[c] # overall kappa
   
   return ddia,ddialo,ddiahi,tkappa,dmass


