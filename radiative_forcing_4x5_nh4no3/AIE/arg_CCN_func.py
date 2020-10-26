from pylab import *
#from numba import autojit

#@autojit
def arg_CCN(temp,pres,V,dN,scrit,scritl):

   # temp in K
   # pres in Pa
   # V in m/s
   # dN in bin in # m-3
   # critical supersaturation for bin mid
   # critical supersaturation for bin limits (nbins+1)

   tempc=temp-273.15
   g = 9.8 # m/s
   R = 8.314 # gas constant J/mol/K
   MWwater = 0.018 # mol wgt water [kg/mol]
   MWair = 0.029 # mol wgt air [kg/mol]
   STwater = 0.072 # surface tension of water [N m-1 or Kg s-2]
   denswater = 1000. # density of water [kg m-3]
   Afactor=4.*STwater*MWwater/(R*temp*denswater) # m
   Lw = (2500.8-2.36*tempc+0.0016*tempc**2-0.00006*tempc**3)*1000. # Latent heat water J/kg
   Cpa=1.01E3 # Heat capacity of air J/kg/K
   Psat=610.78*exp(tempc/(tempc+238.3)*17.2694) # sat vap pres water, Pa
   wvd=-2.775E-6+4.479e-8*temp+1.656e-10*temp**2 # water vap diff [m2/s]
   cona=0.024 # conductivity of air [W/m/K]
   Ntot=dN.sum()
   
   alpha=g*MWwater*Lw/(Cpa*R*temp**2) - g*MWair/(R*temp) # eq 11 in ARG1998
   gamma=R*temp/(Psat*MWwater) + MWwater*Lw**2/(Cpa*pres*MWair*temp) # eq 12 in ARG1998
   #G=4./(denswater*R*temp/(Psat*wvd*MWwater)+Lw*denswater/(cona*temp*(Lw*MWwater/(R*temp)-1.)))
   G=1./(denswater*R*temp/(Psat*wvd*MWwater)+Lw*denswater/(cona*temp)*(Lw*MWwater/(R*temp)-1.))
   # ARG1998 EQ16
   
   zeta=2.*Afactor/3.*(alpha*V/G)**(1./2.) # ARG2002 Eq5
   eta=(alpha*V/G)**(3./2.)/(2.*pi*denswater*gamma*Ntot) # ARG2002 Eq6
   
   Se=(Ntot/(dN/(scrit**(2./3.))).sum())**(3./2.) # ARG2002 Eq8
   Smax=Se/((0.5*(zeta/eta)**(3./2.)+(Se**2/(eta+3.*zeta))**(3./4.))**(1./2.))
   
   ind=where(scritl[:-1]<Smax)[0]
#   print ind
#   print scritl[:-1]
#   print temp, pres, V
   
   if len(ind)>0:
      CDNC=dN[ind].sum() # sum fully activated bins
      if ind[0]>0:
         tie=(Smax-scritl[ind[0]])/(scritl[ind[0]-1]-scritl[ind[0]])
         CDNC=CDNC+tie*dN[ind[0]-1]
   else:
      CDNC = 0.
      
   return Smax,CDNC*1E-6

