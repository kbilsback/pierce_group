 # READ rrtmg_main out file

from pylab import *
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap, addcyclic

#### -------------------------------------------------------------------------------####
#### Routines
#### -------------------------------------------------------------------------------####
def make_map(lat, lon, data):
    m = Basemap(projection='cyl', llcrnrlat=-90, urcrnrlat=90\
                , llcrnrlon=-180, urcrnrlon=180, resolution='c')
    m.drawcoastlines()
    m.drawmapboundary()
    m.drawcountries()
    m.drawparallels([-90,-60,-30,0,30,60,90], linewidth=0.25)
    m.drawmeridians([-120,-60,0,60,120], linewidth=0.25)
    datap, lonp = addcyclic(data, lon)
    LONS, LATS = np.meshgrid(lonp, lat)
    x, y = m(LONS, LATS)
    return m, x, y, datap

def getCLR():
    clr=((0.,0.,.7),
         (0.,0.,1.),
         (0.3,0.5,1.),
         (0.45,0.72,1.),
         (0.6,0.9,1.),
         (1.,1.,1.),
         (1.,0.9,0.6),
         (1.,0.72,0.45),
         (1.,0.5,0.3),
         (1.,0.,0.),
         (0.7,0.,0.))
    
    return clr


#### -------------------------------------------------------------------------------####
#### Read in output file
#### -------------------------------------------------------------------------------####
nc_file = Dataset('/pierce-scratch/kbilsback/geos_chem/projects/china_rescoal/radiative_effects/met_4x5/met_albedo_4x5_w1516.nc')
lat = nc_file.variables['lat'][:]
lon = nc_file.variables['lon'][:]

# Output files
nc1 = Dataset('/pierce-scratch/kbilsback/geos_chem/projects/china_rescoal/radiative_effects/DRE/core-shell/RRTMG_output/globe_off_w1516.nc')
nc2 = Dataset('/pierce-scratch/kbilsback/geos_chem/projects/china_rescoal/radiative_effects/DRE/core-shell/RRTMG_output/globe_w1516.nc')

out = 'globe_w1516'

lv=-1

# READ IN SIMULATION 1
SWUFLUX1 = nc1.variables['SWUFLUX'][:]
SWDFLUX1 = nc1.variables['SWDFLUX'][:]
SWUFLUXC1 = nc1.variables['SWUFLUC'][:]
SWDFLUXC1 = nc1.variables['SWDFLUC'][:]

# FOR SIMULATION 2
SWUFLUX2 = nc2.variables['SWUFLUX'][:]
SWDFLUX2 = nc2.variables['SWDFLUX'][:]
SWUFLUXC2 = nc2.variables['SWUFLUC'][:]
SWDFLUXC2 = nc2.variables['SWDFLUC'][:]

#### -------------------------------------------------------------------------------####
#### Subtract TOA and take average
#### -------------------------------------------------------------------------------####
# NET FLUX ##
net1 = SWDFLUX1[:,:,:] - SWUFLUX1[:,:,:]
net2 = SWDFLUX2[:,:,:] - SWUFLUX2[:,:,:]

netTOA = net1[-1,:,:] - net2[-1,:,:]
netTOAavg = netTOA[:,:].mean(1)
netTOAavg = reshape(netTOAavg, (len(lat), len(lon)))

LONS, LATS = meshgrid(lon, lat)
glob_avg = average(netTOAavg, weights=cos(pi/180.*LATS))

# SAVE OUTPUT
np.savez('/pierce-scratch/kbilsback/geos_chem/projects/china_rescoal/radiative_effects/DRE/core-shell/dre_output/'+out+'_DRE_INT.npz', DRE=netTOAavg, lat=lat,lon=lon)


#### -------------------------------------------------------------------------------####
#### Plotting
#### -------------------------------------------------------------------------------####
levels = array([-0.4,-0.3,-0.2,-0.1,-0.05,-0.01,0.01,0.05,0.1,0.2,0.3,0.4])
#positive_levels = np.logspace(np.log10(1), np.log10(56), 6)
#levels = np.zeros(12)
#for i in range(6):
#    levels[5-i] = -positive_levels[i]
#    levels[6+i] = positive_levels[i]

figure()
m,x,y,datap = make_map(lat,lon,netTOAavg)
m.drawcountries()
overlay = contourf(x,y,datap, levels, colors=getCLR(), extend='both')
underclr = (0.00,0.00,0.20)
overclr = (0.40,0.00,0.00)
overlay.cmap.set_under(underclr)
overlay.cmap.set_over(overclr)
cbar = colorbar(orientation='horizontal')
cbar.set_label('[W m$^{-2}$]')
savefig('/pierce-scratch/kbilsback/geos_chem/projects/china_rescoal/radiative_effects/DRE/core-shell/dre_output/'+out+'_DRE_INT.png')
show()
