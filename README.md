# ROMS NetCDF Cleaner
This code cleans ROMS output NetCDF file; keeps only required parameters.
Note that the path to the ROMS output NetCDF file should be set in the specified line of the code.

The following parameters are extracted and saved into the new NetCDF file:<br>
*u_eastward*<br>
*v_northward*<br>
*w_upward*<br>
*salinity*<br>
*temp (temprature)*<br>
*zeta (absolute sea surface height)*<br>

## Prerequisites
NumPy, SciPy.interpolate and netCDF4 installed either by pip or conda
