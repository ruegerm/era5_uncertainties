
from datetime import datetime, timedelta
from netCDF4 import Dataset, MFDataset, num2date, date2num
import matplotlib.pyplot as plt
import numpy as np
from timeit import default_timer as timer
from os.path import join as path_join
from config_wat_plot import era5_data_dir, plot_output_dir, wind_file_name_format, geopotential_file_name
    

# all era5 altitude levels
altitude_levels = [13077.79, 10986.70, 8951.30, 6915.29, 4892.26, 3087.75, 2653.58, 2260.99, 1910.19,
                   1600.04, 1459.58, 1328.43, 1206.21, 1092.54, 987.00, 889.17, 798.62, 714.94, 637.70,
                   566.49, 500.91, 440.58, 385.14, 334.22, 287.51, 244.68, 205.44, 169.50, 136.62,
                   106.54, 79.04, 53.92, 30.96, 10.00]

# Heights at which the wattisham data is evaluated
wat_heights_meaned= [187, 260, 334, 347, 408, 481, 552, 555, 628, 702, 757, 776, 849, 923, 962, 996, 1070, 1144, 1167, 1217, 1291, 1365, 1371, 1438, 1512, 1576, 1585, 1659]

def get_surface_elevation(wind_lat, wind_lon):
    """Determine surface elevation using ERA5 geopotential data file.

    Args:
        wind_lat (list): Latitudes used in the wind data file.
        wind_lon (list): Longitudes used in the wind data file.

    Returns:
        np.ndarray: Array containing the surface elevation in meters above mean sea level.

    """
    # Load the NetCDF file containing the geopotential of Europe.
    nc = Dataset(path_join(era5_data_dir, geopotential_file_name))
    
    # Read the variables from the netCDF file.
    geopot_lat = nc.variables['latitude'][:]
    geopot_lon = nc.variables['longitude'][:]
    
    
    # Check if wind and geopotential data use same grid.
    assert np.array_equal(geopot_lat, wind_lat) and np.array_equal(geopot_lon, wind_lon), \
        "Requested latitudes and/or longitudes do not correspond to those in the NetCDF file."

    geopot_z = nc.variables['z'][0, :, :]
    nc.close()

    surface_elevation = geopot_z/9.81
    print("Minimum and maximum elevation found are respectively {:.1f}m and {:.1f}m, removing those below zero."
          .format(np.amin(surface_elevation), np.amax(surface_elevation)))

    # Get rid of negative elevation values.
    for i, row in enumerate(surface_elevation):
        for j, val in enumerate(row):
            if val < 0.:
                surface_elevation[i, j] = 0.

    return surface_elevation


def read_era5_data(start_year, final_year):
    # Data from era5 netcdf files for wattisham coordinates and year 2010 is read in and returned as a set of lists
    
    
    # Construct the list of input NetCDF files
    netcdf_files = []
    for y in range(start_year, final_year+1):
        for m in range(1, 13):
            netcdf_files.append(path_join(era5_data_dir, wind_file_name_format.format(y, m)))

    # Load the data from the NetCDF files.
    nc = MFDataset(netcdf_files)
    
    
    
    # ~ # Read the variables from the netCDF file.
    lons = nc.variables['longitude'][:]
    lats = nc.variables['latitude'][:]
    levels = nc.variables['level'][:]  # Model level numbers.
    hours = nc.variables['time'][:]  # Hours since 1900-01-01 00:00:0.0, see: print(nc.variables['time']).
    
    return nc, lons, lats, levels, hours


def process_era5_data(start_year, final_year):
    # creates a netcdf file for era5 data similar to the wattisham netcdf file
    nc, lons, lats, levels, hours = read_era5_data(start_year, final_year)
    print('ERA5 Data has been read in')
    
    print('Getting surface elevation...')
    surface_elevation = get_surface_elevation(lats, lons)
    print('Got surface elevation. Now getting speed values...')
    v_levels_east = nc.variables['u'][:,:,:,:]
    v_levels_north = nc.variables['v'][:,:,:,:]
    v_levels = (v_levels_east**2 + v_levels_north**2)**.5
    print('Got speed values.')
    
    #create netcdf file for processed data to be written to
    nc_out = Dataset('era5_data_sorted.nc', "w", format="NETCDF3_64BIT_OFFSET")
    nc_out.createDimension("length", len(hours)*len(wat_heights_meaned))
    
    print('dimension length', len(hours)*len(wat_heights_meaned))
    
    nc_date = nc_out.createVariable("date", "i4", ("length",))
    nc_time = nc_out.createVariable("time", "i4", ("length",))
    nc_hour = nc_out.createVariable("hour", "i4", ("length",))
    # ~ nc_lat = nc_out.createVariable("lat", "f4", ("length",))
    # ~ nc_lon = nc_out.createVariable("lon", "f4", ("length",))
    nc_speed = nc_out.createVariable("speed", "f4", ("length",))
    nc_height = nc_out.createVariable("height", "i4", ("length",))
    nc_wat_height = nc_out.createVariable("wat_height", "i4", ("length",))
    
    # create dummy arrays 
    hourss=[]
    dates=[]
    times=[]
    # ~ latss=[]
    # ~ lonss=[]
    speeds= []
    heights= []
    wat_heights_meaneds= []
    
    
    print('Now starting processing loop...')
    # ~ #loop over the input files variables
    for i_hr,hour in enumerate(hours):
        
        # only select the coordinates closest to Wattisham airfield, lat=52 (row in matrix = 6), lon=1 (row in matrix = 22)
        
        altitudes_of_interest = wat_heights_meaned + surface_elevation[6,22] # calculate heights_over_ground from model_levels at wattisham airfield
        sp_req_alt = np.zeros((len(hours), len(altitudes_of_interest)))  # result array for writing interpolated data
        sp_req_alt[i_hr, :] = np.interp(altitudes_of_interest, altitude_levels[::-1],v_levels[i_hr, ::-1, 6, 22]) # write interpolated hgt_over_grnd into result array
             
        for i_alt, alt in enumerate(altitudes_of_interest):
            # write the actual data of every altitude of every hour into arrays
            hourss.append(hour)
            # ~ latss.append(v_lat)
            # ~ lonss.append(v_lon)
            speeds.append(sp_req_alt[i_hr,i_alt])
            heights.append(alt)
            wat_heights_meaneds.append(wat_heights_meaned[i_alt])
                      
    print('Now converting time:')
    # convert hours-since-1900 to date and time
    datetime= num2date(hourss[:],units='hours since 1900-01-01 00:00:0.0.')
    
    print('converting done')
    dates = [a.strftime("%Y%m%d")[2:] for a in datetime]
    
    times = [int(b.strftime("%H")) for b in datetime]
    
    # write data into nc_variables
    nc_hour[:]=hourss
    nc_date[:]=dates
    nc_time[:]=times
    # ~ nc_lat[:]=latss
    # ~ nc_lon[:]=lonss
    nc_speed[:]=speeds
    nc_height[:]=heights             
    nc_wat_height[:]=wat_heights_meaneds
    
    
    nc_out.close()       
         
    print('nc file done')
        
    
def main():
    
     
    process_era5_data(2010,2010)
    
    
    
    
           
if __name__ == "__main__":
    main()
