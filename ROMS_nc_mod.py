import numpy as np
import scipy.interpolate
import netCDF4
import datetime as dt

strt_time = dt.datetime.now()

"""
Initializing Variables
"""
lon_rho = []
lat_rho = []
s_rho = []
times = []

"""
Reading Source File Data
"""
dataset = netCDF4.Dataset('his_00002.nc', 'r')      # ROMS Output NetCDF File
lon_rho = np.array(dataset.variables['lon_rho'][:])
lat_rho = np.array(dataset.variables['lat_rho'][:])
s_rho[:] = dataset.variables['s_rho'][:]
times[:] = dataset.variables['ocean_time'][:]
times_units = dataset.variables['ocean_time'].units

"""
Creating Destination Grid
"""
lons = np.arange(46.03, 61.15+.03, 0.03)
lats = np.arange(20.91, 30.93+.03, 0.03)
lon_grid,lat_grid = np.meshgrid(lons, lats)

"""
Creating NetCDF File
"""
# Creating NetCDF file
test_file = netCDF4.Dataset('test.nc','w', format='NETCDF4')

# Creating dimensions
test_file.createDimension('lon', len(lons))
test_file.createDimension('lat', len(lats))
test_file.createDimension('z', len(s_rho))
test_file.createDimension('time', None)

# Creating global attributes
test_file.title='Trimmed and Regridded Output File'
test_file.source='ROMS'
test_file.author='Farrokh A. Ghavanini'
test_file.contact='ghavanini[at]gmail[dot]com'

# Creating variables
lon = test_file.createVariable('lon', np.float64, ('lon',))
lon.units = 'degrees_east'
lon.long_name = 'longitude'
lon.point_spacing = 'even'
lat = test_file.createVariable('lat', np.float64, ('lat',))
lat.units = 'degrees_north'
lat.long_name = 'latitude'
lat.point_spacing = 'even'
z = test_file.createVariable('z', np.float64, ('z',))
z.long_name = 'Z-Level'
z.point_spacing = 'even'
time = test_file.createVariable('time', np.float64, ('time',))
time.units = times_units
time.calendar = 'gregorian'
time.long_name = 'time'
u_eastward = test_file.createVariable('u_eastward', np.float64,('time','z','lat','lon'), zlib=True)
u_eastward.units = 'm s-1'
u_eastward.long_name = 'eastward_sea_water_velocity'
v_northward = test_file.createVariable('v_northward', np.float64,('time','z','lat','lon'), zlib=True)
v_northward.units = 'm s-1'
v_northward.long_name = 'northward_sea_water_velocity'
w_upward = test_file.createVariable('w', np.float64,('time','z','lat','lon'), zlib=True)
w_upward.units = 'm s-1'
w_upward.long_name = 'upward_sea_water_velocity'
salinity = test_file.createVariable('salt', np.float64,('time','z','lat','lon'), zlib=True)
salinity.long_name = 'salinity'
temp = test_file.createVariable('temp', np.float64,('time','z','lat','lon'), zlib=True)
temp.units = 'Celsius'
temp.long_name = 'potential temperature'
zeta = test_file.createVariable('zeta', np.float64,('time','lat','lon'), zlib=True)
zeta.units = 'meters'
zeta.long_name = 'free-surface'

# Filling variables
lon[:] = lons
lat[:] = lats
z[:] = s_rho
time[:] = times
u_eastward[:,:,:,:] = np.empty((0, len(s_rho), len(lats), len(lons)))
v_northward[:,:,:,:] = np.empty((0, len(s_rho), len(lats), len(lons)))
w_upward[:,:,:,:] = np.empty((0, len(s_rho), len(lats), len(lons)))
salinity[:,:,:,:] = np.empty((0, len(s_rho), len(lats), len(lons)))
temp[:,:,:,:] = np.empty((0, len(s_rho), len(lats), len(lons)))
zeta[:,:,:] = np.empty((0, len(lats), len(lons)))
lon_array = lon_rho.flatten()
lat_array = lat_rho.flatten()
for timestep, time_value in enumerate(times):
    print('\n-> Processing Timestep No. ' + str(timestep))
    z_array = (np.array(dataset.variables['zeta'][timestep, :, :])).flatten()
    z = scipy.interpolate.griddata((lon_array,lat_array), z_array, (lon_grid,lat_grid), method='linear')
    z[z>5] = np.NaN
    z[z<-5] = np.NaN
    zeta[timestep, :,:] = z
    for z_level, z_value in enumerate(s_rho):
        u_array = (np.array(dataset.variables['u_eastward'][timestep, z_level, :, :])).flatten()
        u = scipy.interpolate.griddata((lon_array,lat_array), u_array, (lon_grid,lat_grid), method='linear')
        u[u>10] = np.NaN
        u[u<-10] = np.NaN
        u_eastward[timestep, z_level, :,:] = u
        v_array = (np.array(dataset.variables['v_northward'][timestep, z_level, :, :])).flatten()
        v = scipy.interpolate.griddata((lon_array,lat_array), v_array, (lon_grid,lat_grid), method='linear')
        v[v>10] = np.NaN
        v[v<-10] = np.NaN
        v_northward[timestep, z_level, :,:] = v        
        w_array = (np.array(dataset.variables['w'][timestep, z_level, :, :])).flatten()
        w = scipy.interpolate.griddata((lon_array,lat_array), w_array, (lon_grid,lat_grid), method='linear')
        w[w>10] = np.NaN
        w[w<-10] = np.NaN
        w_upward[timestep, z_level, :,:] = w
        s_array = (np.array(dataset.variables['salt'][timestep, z_level, :, :])).flatten()
        s = scipy.interpolate.griddata((lon_array,lat_array), s_array, (lon_grid,lat_grid), method='linear')
        s[s>250] = np.NaN
        salinity[timestep, z_level, :,:] = s
        t_array = (np.array(dataset.variables['temp'][timestep, z_level, :, :])).flatten()
        t = scipy.interpolate.griddata((lon_array,lat_array), t_array, (lon_grid,lat_grid), method='linear')
        t[t>50] = np.NaN
        t[t<-20] = np.NaN
        temp[timestep, z_level, :,:] = t
        print(' --> Level ' + str(z_level) + ' Processed.')
    print('\n --> Timestep No. ' + str(timestep) + ' Processed. (ET: ' + str((dt.datetime.now() - strt_time).seconds) + 's)')    
    if timestep == 1:
        break

# Closing NetCDF file
test_file.close()