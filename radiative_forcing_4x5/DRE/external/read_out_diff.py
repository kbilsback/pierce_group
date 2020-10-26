# Plot RRTMG output
# Kelsey Bilsback January 2018

#Hack to fix missing PROJ4 env var
import os
import conda

conda_file_dir = conda.__file__
conda_dir = conda_file_dir.split('lib')[0]
proj_lib = os.path.join(os.path.join(conda_dir, 'share'), 'proj')
os.environ["PROJ_LIB"] = proj_lib


from pylab import *
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import BoundaryNorm

#### -------------------------------------------------------------------------------####
#### Routines
#### -------------------------------------------------------------------------------####
def make_map(lat, lon):
    m = Basemap(projection = 'cyl',
                llcrnrlat = 29.0,
                urcrnrlat = 43.75,
                llcrnrlon = 108.75,
                urcrnrlon = 131.875,
                resolution = 'i')
    
    m.drawcoastlines()
    m.drawmapboundary()
    m.drawcountries()
    #m.drawparallels([-90,-60,-30,0,30,60,90], linewidth=0.25)
    #m.drawmeridians([-120,-60,0,60,120], linewidth=0.25)
    m.readshapefile('/pierce-scratch/kbilsback/geos_chem/projects/china_rescoal/shapefiles/gadm36_CHN_1', 'prov', linewidth=0.25)
    #m.readshapefile('/pierce-scratch/kbilsback/geos_chem/projects/china_rescoal/shapefiles/gadm36_CHN_2', 'dist', linewidth=0.5)
    #m.readshapefile('/pierce-scratch/kbilsback/geos_chem/projects/china_rescoal/shapefiles/gadm36_CHN_3', 'town', linewidth=0.15)
    
    LONS, LATS = np.meshgrid(lon, lat)
    x, y = m(LONS, LATS)
    return m, x, y

def get_cmap():
    
    colors=[(0.,0.,.7),
            (0.,0.,1.),
            (0.3,0.5,1.),
            (0.45,0.72,1.),
            (0.6,0.9,1.),
            (1.,1.,1.),
            (1.,0.9,0.6),
            (1.,0.72,0.45),
            (1.,0.5,0.3),
            (1.,0.,0.),
            (0.7,0.,0.)]

    c_map = matplotlib.colors.LinearSegmentedColormap.from_list("", colors)
    
    return c_map


#### -------------------------------------------------------------------------------####
#### Read in output file
#### -------------------------------------------------------------------------------####
nc_file = Dataset('/pierce-scratch/kbilsback/geos_chem/projects/china_rescoal/radiative_effects_025x03125/met_025x03125/met_albedo_025x03125_w1415.nc')
lat = nc_file.variables['lat'][:]
lon = nc_file.variables['lon'][:]

# Output files
nc1 = Dataset('/pierce-scratch/kbilsback/geos_chem/projects/china_rescoal/radiative_effects_025x03125/DRE/external/RRTMG_output/china_boff_w1415.nc')
nc2 = Dataset('/pierce-scratch/kbilsback/geos_chem/projects/china_rescoal/radiative_effects_025x03125/DRE/external/RRTMG_output/china_w1415.nc')

out = 'china_bthoff'

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
np.savez('/pierce-scratch/kbilsback/geos_chem/projects/china_rescoal/radiative_effects_025x03125/AIE/aie_output/'+out+'_AIE_w1415.npz', AIE=netTOAavg, lat=lat,lon=lon)


#### -------------------------------------------------------------------------------####
#### Plotting
#### -------------------------------------------------------------------------------####

# define color levels
levels = array([-0.4,-0.3,-0.2,-0.1,-0.05,-0.01,0.01,0.05,0.1,0.2,0.3,0.4])

close('all')
# make figure
fig = figure(figsize=(8, 6), dpi=120, facecolor='w', edgecolor='w')

m,x,y = make_map(lat,lon)
c_map = get_cmap()
norm = BoundaryNorm(levels, ncolors=c_map.N, clip=True)

p = m.pcolormesh(x, y, netTOAavg, cmap=c_map, norm=norm)

# astetics 
cbar = colorbar(orientation='horizontal', ticks=levels, extend='both')
cbar.set_label('[$\Delta$W m$^{-2}$]')
cbar.ax.set_xticklabels(['-0.4','-0.3','-0.2','-0.1','-0.05','-0.01','0.01','0.05','0.1','0.2','0.3','0.4']) 

# save and close
savefig('/pierce-scratch/kbilsback/geos_chem/projects/china_rescoal/radiative_effects_025x03125/AIE/aie_output/'+out+'_AIE_w1415.png')
show()

close('all')
