from netCDF4 import Dataset, MFDataset
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats 
from timeit import default_timer as timer
from os.path import join as path_join
from config_wat_plot import era5_data_dir, plot_output_dir, wind_file_name_format
from config_wat_analyse import output_nc_file
from wat_analyse import sort_by_index, read_out_data, get_distance_directions    


req_altitudes = [184, 256, 328, 401,473, 545, 617, 690, 762, 834, 907, 979, 
                       1051, 1124, 1196, 1268, 1341, 1413, 1485, 1557, 1630, 1702]

heights_lidar1_mean = [187, 260, 334, 408, 481, 555, 628, 702, 776, 849, 923, 996, 1070, 1144, 1217, 1291, 1365, 1438, 1512, 1585, 1659]
heights_lidar2 = [347, 552, 757, 962, 1167, 1371, 1576]
heights_meaned= [187, 260, 334, 347, 408, 481, 552, 555, 628, 702, 757, 776, 849, 923, 962, 996, 1070, 1144, 1167, 1217, 1291, 1365, 1371, 1438, 1512, 1576, 1585, 1659]
times=[i for i in range(24)]


def read_wat_data(nc_input_wat):
    # Data from the Wattisham netcdf file is read in and returned as a set of lists
    
    nc = Dataset(nc_input_wat,'r')
    
    length= len(nc.dimensions['length'])
        
    dates = nc.variables["date"][:]
    times = nc.variables["time"][:]
    times_seconds = nc.variables["time_seconds"][:]
    heights = nc.variables["height"][:]
    speeds_ur = nc.variables["speed"][:]
    directions = nc.variables["direction"][:]
    vel_2ds_ur = nc.variables["vel_2d"][:]
    vel_us_ur = nc.variables["vel_u"][:]
    vel_vs_ur = nc.variables["vel_v"][:]
    vel_ws_ur = nc.variables["vel_w"][:]
    met_qcs = nc.variables["met_qc"][:]
    heights_mean = nc.variables["height_mean"][:]
    
    # round the values in the list
    
    speeds = [ round(elem, 3) for elem in speeds_ur]
    vel_2ds = [ round(elem, 3) for elem in vel_2ds_ur]
    vel_us = [ round(elem, 3) for elem in vel_us_ur]
    vel_vs = [ round(elem, 3) for elem in vel_vs_ur]
    vel_ws = [ round(elem, 3) for elem in vel_ws_ur]
    
    all_lines = nc_input_to_list(length, dates, times, times_seconds, heights, speeds, directions, vel_us, vel_vs, vel_ws, met_qcs, heights_mean, vel_2ds)
    
    
    print('Wattisham Data is read in')
    return all_lines

def read_wat_era5_data(input_file):
    # Data from the processed era5 netcdf file is read in and returned as a set of lists
    
    nc = Dataset(input_file,'r')
    
    length= len(nc.dimensions['length'])
    
    hours = nc.variables["hour"][:]    
    dates = nc.variables["date"][:]
    times = nc.variables["time"][:]
    # ~ lats = nc.variables["lat"][:]
    # ~ lons = nc.variables["lon"][:]
    heights = nc.variables["height"][:]
    wat_heights = nc.variables["wat_height"][:]
    speeds = nc.variables["speed"][:]
    
    # round the values in the lists
    
    # ~ lats = [ round(elem, 3) for elem in lats]
    # ~ lons = [ round(elem, 3) for elem in lons]
    speeds = [ round(elem, 2) for elem in speeds]
    
    all_lines=[]
    dum=[]
    for a in range(length):
        dum.append(hours[a])
        dum.append(dates[a])
        dum.append(times[a])
        # ~ dum.append(lats[a])
        # ~ dum.append(lons[a])
        dum.append(heights[a])
        dum.append(wat_heights[a])
        dum.append(speeds[a])
        
        all_lines.append(dum)
        dum=[]
        
    return all_lines
    print('Era5 Data is read in')



def plot_hist(i_figure, input_hist, savefile, origin_hist_1=None, name_1='input_1',xlabel=None, ylabel='events',nbins=None,binrange=None,title='',histtype='stepfilled'):
    # function for plotting one hist on the canvas
    # name, axis labels, nbins, binrange, title and histtype can be specified
    
    plt.figure(i_figure)
    plt.hist(input_hist,nbins,binrange,histtype='stepfilled')
    plt.xlabel(xlabel, fontsize=18)
    plt.ylabel(ylabel, fontsize=18)
    
    #calculate mean and stdev
    a=np.array(input_hist)
    mean=a.mean()
    stdev=np.std(a)
    rms = np.sqrt(np.mean(np.square(a)))
    
    error_rel1=None
    
    if origin_hist_1!= None:
        b=np.array(origin_hist_1)
        bmean=b.mean()
        error_rel1 = round(rms/bmean,2)
    
    plt.annotate('{}'.format(name_1), xy=(0.8, 0.95), xycoords='axes fraction',color='b')
    plt.annotate('Events: {}'.format(len(a)), xy=(0.8, 0.9), xycoords='axes fraction',color='b')
    plt.annotate('Mean: {:.2f}'.format(mean), xy=(0.8, 0.85), xycoords='axes fraction',color='b')
    plt.annotate('RMS: {:.2f}'.format(rms), xy=(0.8, 0.8), xycoords='axes fraction',color='b')
    if error_rel1 != None:
        plt.annotate('RMS/mean velocity: {:.2f}'.format(error_rel), xy=(0.8, 0.75), xycoords='axes fraction',color='b')
    # ~ plt.text(30, max(a), 'Mean: {:.2f}'.format(mean))
    plt.title(title, fontsize = 18)
    plt.savefig(savefile)
    plt.close()

def plot_hist_2(i_figure, input_hist_1, input_hist_2, savefile, origin_hist_1=None, name_1='input_1', name_2='input_2',style='',xlabel=None, ylabel='events',nbins=100,binrange=None,title='',histtype='step'):
    #function for plotting two hists onto same canvas
    
    
    #calculate mean and rms
    a=np.array(input_hist_1)
    mean_1=a.mean()
    rms_1 = np.sqrt(np.mean(np.square(a)))
    b=np.array(input_hist_2)
    mean_2=b.mean()
    rms_2 = np.sqrt(np.mean(np.square(b)))
    error_rel1= round(rms_1/mean_1,2)
    error_rel2= round(rms_1/mean_2,2)
    
    plt.figure(i_figure)
    plt.hist(input_hist_1,nbins,binrange,histtype='step',color='b')
    plt.hist(input_hist_2,nbins,binrange,histtype='step',color='g')

    plt.xlabel(xlabel, fontsize= 18)
    plt.ylabel(ylabel, fontsize= 18)
    
    error_rel1=None
    
    if origin_hist_1!= None:
        c=np.array(origin_hist_1)
        bmean=c.mean()
        error_rel1 = round(rms_1/cmean,2)
    
    
    
    plt.annotate('{}'.format(name_1), xy=(0.8, 0.95), xycoords='axes fraction', color='b')
    plt.annotate('Events: {}'.format(len(input_hist_1)), xy=(0.8, 0.9), xycoords='axes fraction', color='b')
    plt.annotate('Mean: {:.2f}'.format(mean_1), xy=(0.8, 0.85), xycoords='axes fraction', color='b')
    plt.annotate('RMS: {:.2f}'.format(rms_1), xy=(0.8, 0.8), xycoords='axes fraction', color='b')
    if error_rel1 != None:
        plt.annotate('RMS/mean velocity: {:.2f}'.format(error_rel), xy=(0.8, 0.75), xycoords='axes fraction',color='b')
    plt.annotate('{}'.format(name_2), xy=(0.8, 0.7), xycoords='axes fraction', color='g')
    plt.annotate('Events: {}'.format(len(input_hist_2)), xy=(0.8, 0.65), xycoords='axes fraction', color='g')
    plt.annotate('Mean: {:.2f}'.format(mean_2), xy=(0.8, 0.6), xycoords='axes fraction', color='g')
    plt.annotate('RMS: {:.2f}'.format(rms_2), xy=(0.8, 0.55), xycoords='axes fraction', color='g')
    # ~ plt.text(30, max(a), 'Mean: {:.2f}'.format(mean))
    plt.title(title, fontsize=18)
    plt.savefig(savefile)
    plt.close()
      
def plot_scatter_2D(i_figure, input_hist_1, input_hist_2, savefile, name_1, name_2,style='',xlabel=None, ylabel=None,title=''):
    
    plt.figure(i_figure)
    # Calculate the point density
    xy = np.vstack([input_hist_1,input_hist_2])
    z = stats.gaussian_kde(xy)(xy)
    plt.scatter(input_hist_1,input_hist_2,s=10,c=z)

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    
    # ~ plt.text(30, max(a), 'Mean: {:.2f}'.format(mean))
    plt.title(title)
    plt.savefig(savefile)
    plt.close()


def plot_delta_dir_time(input_all_lines,heights, name_input_1='input_1', output_path='./plots/testplots', index=5,Nbins=100,Binrange=(-2,35)):
    
    hist_dummy=[]
    for i_time,time in enumerate(times):
        if i_time==0:
            continue
        for c,v in enumerate(input_all_lines):
            
            if v[1] == time and input_all_lines[c-1][1]== times[i_time-1] :
                
                hist_dummy.append( get_distance_directions(v[index],input_all_lines[c-1][index]) )
        
        output= '{}/delta_winddirection_time_{}_to_{}'.format(output_path, time,times[i_time-1])
        titel='time_{}_to_{}'.format(time,times[i_time-1],len(hist_dummy)) 
        plot_hist(i_time, hist_dummy, savefile=output, name_1=name_input_1, xlabel='dir2-dir1',nbins=Nbins,binrange=Binrange,title=titel )
        print('plot {} done'.format(i_time))        
        hist_dummy=[]

def plot_delta_dir_time_2(input_all_lines_1,input_all_lines_2,heights, name_input_1='input_1', name_input_2='input_2',output_path='./plots/testplots',Nbins=100,Binrange=(-2,35),index=5):
    
    hist_dummy_1=[]
    hist_dummy_2=[]
    
    for i_time,time in enumerate(times):
       
        if i_time==0:
            continue
        for c,v in enumerate(input_all_lines_1):
            
            if v[1] == time and input_all_lines_1[c-1][1]== times[i_time-1] :
                
                hist_dummy_1.append( get_distance_directions(v[index],input_all_lines_1[c-1][index]) )
        
        for c,v in enumerate(input_all_lines_2):
            
            if v[1] == time and input_all_lines_2[c-1][1]== times[i_time-1] :
                
                hist_dummy_2.append( get_distance_directions(v[index],input_all_lines_2[c-1][index]) )
        
        output= '{}/delta_winddirection_time_{}_to_{}'.format(output_path,time,times[i_time-1])
        titel='time_{}_to_{}'.format(time,times[i_time-1]) 
        plot_hist_2(i_time, hist_dummy_1,hist_dummy_2, savefile=output, name_1=name_input_1, name_2=name_input_2, xlabel='dir2-dir1',nbins=Nbins,binrange=Binrange,title=titel )
        print('plot {} done'.format(i_time))        
        hist_dummy_1=[]
        hist_dummy_2=[]
    
                    
def plot_delta_dir_height(input_all_lines, heights, origin_hist1=None, name_input_1='input_1', output_path='./plots/testplots',index=5,Nbins=100,Binrange=(-2,35)):
    
    hist_dummy=[]
    for i_height,height in enumerate(heights):
        
        if i_height==0:
            continue
            
        for c,v in enumerate(input_all_lines):
            
            if input_all_lines[c][10] == height and ( input_all_lines[c-1][10]== heights[i_height-1] ):
                
                hist_dummy.append( get_distance_directions(v[index],input_all_lines[c-1][index]) )
                
        output= '{}/delta_winddirection_height_{}_to_{}'.format(output_path,height,heights[i_height-1])
        titel='height_{}_to_{}'.format(height,heights[i_height-1],len(hist_dummy)) 
        plot_hist(i_height, hist_dummy, origin_hist_1= origin_hist1, savefile=output, name_1=name_input_1, xlabel='dir2-dir1',nbins=Nbins,binrange=Binrange,title=titel )
        print('plot {} done'.format(i_height))        
        hist_dummy=[]  
        
def plot_delta_dir_height_2(input_all_lines_1, input_all_lines_2,heights, name_input_1='input_1', name_input_2='input_2', output_path='./plots/testplots',Nbins=100,Binrange=(-2,35),index=5):
    
    hist_dummy_1=[]
    hist_dummy_2=[]
    
    
    for i_height,height in enumerate(heights):
        
        if i_height==0:
            continue
            
        for c,v in enumerate(input_all_lines_1):
            
            if input_all_lines_1[c][10] == height and ( input_all_lines_1[c-1][10]== heights[i_height-1] ):
                
                hist_dummy_1.append( get_distance_directions(v[index],input_all_lines_1[c-1][index]) )
        
        for c,v in enumerate(input_all_lines_2):
            
            if input_all_lines_2[c][10] == height and ( input_all_lines_2[c-1][10]== heights[i_height-1] ):
                
                hist_dummy_2.append( get_distance_directions(v[index],input_all_lines_2[c-1][index]) )
                
        output= '{}/delta_winddirection_height_{}_to_{}'.format(output_path,height,heights[i_height-1])
        titel='Height: {}m to {}m'.format(height,heights[i_height-1]) 
        plot_hist_2(i_height, hist_dummy_1, hist_dummy_2, savefile=output, name_1=name_input_1, name_2=name_input_2, xlabel='dir2 - dir1',nbins=Nbins,binrange=Binrange,title=titel )
        print('plot {} done'.format(i_height))        
        hist_dummy_1=[]
        hist_dummy_2=[] 
        
        
def plot_delta_vel_time(input_all_lines,times, name_input_1='input_1', output_path='./plots/testplots', index=11,Nbins=100,Binrange=(-2,35)):
    
    hist_dummy=[]
    for i_time,time in enumerate(times):
        if i_time==0:
            continue
        for c,v in enumerate(input_all_lines):
            
            if input_all_lines[c][1] == time and input_all_lines[c-1][1]== times[i_time-1]:
                hist_dummy.append( v[index]-input_all_lines[c-1][index] )
        
        output= '{}/delta_windspeed_time_{}_to_{}'.format(output_path,time,times[i_time-1])
        titel='time_{}_to_{}'.format(time,times[i_time-1],len(hist_dummy)) 
        plot_hist(i_time, hist_dummy, savefile=output, name_1=name_input_1, xlabel='v2-v1',nbins=Nbins,binrange=Binrange,title=titel )
        print('plot {} done'.format(i_time))        
        hist_dummy=[]
        
def plot_delta_vel_time_2(input_all_lines_1, input_all_lines_2,times, name_input_1='input_1', name_input_2='input_2',output_path='./plots/testplots',Nbins=100,Binrange=(-2,35),index=11):
    
    hist_dummy_1=[]
    hist_dummy_2=[]
    for i_time,time in enumerate(times):
        if i_time==0:
            continue
        for c,v in enumerate(input_all_lines_1):
            
            if input_all_lines_1[c][1] == time and input_all_lines_1[c-1][1]== times[i_time-1]:
                hist_dummy_1.append( v[index]-input_all_lines_1[c-1][index] )
        
        for c,v in enumerate(input_all_lines_2):
            
            if input_all_lines_2[c][1] == time and input_all_lines_2[c-1][1]== times[i_time-1]:
                hist_dummy_2.append( v[index]-input_all_lines_2[c-1][index] )
        
        output= '{}/delta_windspeed_time_{}_to_{}'.format(output_path,time,times[i_time-1])
        titel='time_{}_to_{}'.format(time,times[i_time-1]) 
        plot_hist_2(i_time, hist_dummy_1,hist_dummy_2, savefile=output, name_1=name_input_1, name_2=name_input_2, xlabel='v2-v1',nbins=Nbins,binrange=Binrange,title=titel )
        print('plot {} done'.format(i_time))        
        hist_dummy_1=[]
        hist_dummy_2=[]

def plot_delta_vel_height(input_all_lines,heights, origin_hist1=None, name_input_1='input_1', output_path='./plots/testplots', index=11,Nbins=100,Binrange=(-2,35)):
    
    hist_dummy=[]
    for i_height,height in enumerate(heights):
        
        if i_height==0:
            continue
        
        for c,v in enumerate(input_all_lines):
            
            if input_all_lines[c][10] == height and input_all_lines[c-1][10]== heights[i_height-1]:
                
                if v[index]-input_all_lines[c-1][index] > 4:
                    
                    print('{}  {}  {}  '.format(v[0], v[1], v[10] ) )
                    print('{}  {}  {}  difference: {}'.format(input_all_lines[c-1][0], input_all_lines[c-1][1], input_all_lines[c-1][10], v[index]-input_all_lines[c-1][index] ))
                    print('---------------------------------------------')
                    
                hist_dummy.append( v[index]-input_all_lines[c-1][index] )
        
        
        output= '{}/delta_windspeed_height_{}_to_{}'.format(output_path,height,heights[i_height-1])
        titel='Height: {}m to {}m'.format(height,heights[i_height-1],len(hist_dummy)) 
        plot_hist(i_height, hist_dummy, origin_hist_1=origin_hist1,savefile=output, name_1=name_input_1, xlabel='v2-v1',nbins=Nbins,binrange=Binrange,title=titel)
        print('plot {} done'.format(i_height))        
        hist_dummy=[]  

def plot_delta_vel_height_2(input_all_lines_1, input_all_lines_2, heights, name_input_1='input_1', name_input_2='input_2',output_path='./plots/testplots',Nbins=100,Binrange=(-2,35), index=11):
    
    hist_dummy_1=[]
    hist_dummy_2=[]
    for i_height,height in enumerate(heights):
        
        if i_height==0:
            continue
        
        for c,v in enumerate(input_all_lines_1):
            
            if input_all_lines_1[c][10] == height and input_all_lines_1[c-1][10]== heights[i_height-1]:
                
                hist_dummy_1.append( v[index]-input_all_lines_1[c-1][index] )
        
        for c,v in enumerate(input_all_lines_2):
            
            if input_all_lines_2[c][10] == height and input_all_lines_2[c-1][10]== heights[i_height-1]:
                
                hist_dummy_2.append( v[index]-input_all_lines_2[c-1][index] )
        
    
        output= '{}/delta_windspeed_height_{}_to_{}'.format(output_path, height,heights[i_height-1])
        titel='Height: {}m to {}m'.format(height,heights[i_height-1]) 
        plot_hist_2(i_height, hist_dummy_1,hist_dummy_2, savefile=output, name_1=name_input_1, name_2=name_input_2, xlabel='v2 - v1'.format(height,heights[i_height-1]),nbins=Nbins,binrange=Binrange,title=titel)
        print('plot {} done'.format(i_height))        
        hist_dummy_1=[]
        hist_dummy_2=[]
        
def plot_watera5_delta_vel_height(input_wat, input_era5, heights, name_input_1='input_1', output_path='./plots/testplots',Nbins=100,Binrange=(-2,35)):
    
    hist_dummy=[]
    for i_height,height in enumerate(heights):
        
        for c,v in enumerate(input_wat):
            if v[10] == height:
                
                for a,b in enumerate( input_era5[i_height*8759:(i_height*8759+8759)] ):
                    
                    if v[0] == b[1] and v[1]==b[2] and b[4] == height:
                        hist_dummy.append( v[11] - b[5] )  
                    else:
                        continue
                    
            else:
                continue
        
           
        output= '{}/delta_speed_height_wat_to_era5_height{}'.format(output_path,height)
        titel='Height: {}m'.format(height) 
        plot_hist(i_height, hist_dummy, savefile=output, name_1=name_input_1, xlabel='v_wat - v_ERA5',nbins=Nbins,binrange=Binrange,title=titel)
        print('plot {} done'.format(i_height))        
        hist_dummy=[]

      
def plot_watera5_deltarel_vel_height(input_wat, input_era5, heights, name_input_1='input_1', output_path='./plots/testplots',Nbins=100,Binrange=(-2,35)):
    
    hist_dummy=[]
    for i_height,height in enumerate(heights):
        
        for c,v in enumerate(input_wat):
            if v[10] == height:
                
                for a,b in enumerate( input_era5[i_height*8759:(i_height*8759+8759)] ):
                    
                    if v[0] == b[1] and v[1]==b[2] and b[4] == height:
                        hist_dummy.append( (v[11] - b[5] ) / v[11] )  
                    else:
                        continue
                    
            else:
                continue
        
           
        output= '{}/deltarel_speed_height_wat_to_era5_height{}'.format(output_path,height)
        titel='height_{}'.format(height) 
        plot_hist(i_height, hist_dummy, savefile=output, name_1=name_input_1, xlabel='(v_wat-v_ERA5)/v_wat',nbins=Nbins,binrange=Binrange,title=titel)
        print('plot {} done'.format(i_height))        
        hist_dummy=[]        

def plot_watera5_vel_height_2(input_wat, input_era5, heights, name_input_1='input_1', name_input_2='input_2', output_path='./plots/testplots',Nbins=100,Binrange=(-2,35)):
    #   Input : wat_hist, era5_hist, heights_list; Optional: output, Nbins, Binrange
    hist_wat=[]
    hist_era=[]
    
    #   loop over all the heights, each time fill hists with the data and plot the result
    for i_height,height in enumerate(heights):
        
        for c,v in enumerate(input_wat):
            if v[10] == height:
                hist_wat.append(v[11])
        for c,v in enumerate(input_era5):
            if v[4] == height:
                hist_era.append(v[5])
            
        
           
        output= '{}/speed_height_wat_era5_height{}'.format(output_path,height)
        titel='Height: {}m'.format(height) 
        plot_hist_2(i_height, hist_wat, hist_era, savefile=output, name_1=name_input_1, name_2=name_input_2, xlabel='windspeed',nbins=Nbins,binrange=Binrange,title=titel,)
        print('plot {} done'.format(i_height))        
        hist_wat=[] 
        hist_era=[]

def plot_watera5_vel_height_2_same_times(input_wat, input_era5, heights,name_input_1='input_1', name_input_2='input_2',output_path='./plots/testplots',Nbins=100,Binrange=(-2,35)):
    # i plot only those points where there is data in era5 AND wat
    
    hist_wat=[]
    hist_era=[]
    
    for i_height,height in enumerate(heights):
        
        for c,v in enumerate(input_wat):
            if v[10] == height:
                
                for a,b in enumerate( input_era5[i_height*8759:(i_height*8759+8759)] ):
                    
                    if v[0] == b[1] and v[1]==b[2] and b[4] == height:
                        hist_wat.append(v[11])
                        hist_era.append(b[5])  
                    else:
                        continue
                    
            else:
                continue
        
          
        output= '{}/speed_height_wat_era5_sametimes_height{}'.format(output_path,height)
        titel='Height: {}m'.format(height) 
        plot_hist_2(i_height, hist_wat, hist_era, savefile=output, name_1=name_input_1, name_2=name_input_2, xlabel='windspeed [m/s]',nbins=Nbins,binrange=Binrange,title=titel,)
        print('plot {} done'.format(i_height))        
        hist_wat=[] 
        hist_era=[]  
           
def plot_watera5_vel_height_2_same_times_2D(input_wat, input_era5, heights,name_input_1='input_1', name_input_2='input_2',output_path='./plots/testplots'):
    
    hist_wat=[]
    hist_era=[]
    
    for i_height,height in enumerate(heights):
        
        for c,v in enumerate(input_wat):
            if v[10] == height:
                
                for a,b in enumerate( input_era5[i_height*8759:(i_height*8759+8759)] ):
                    
                    if v[0] == b[1] and v[1]==b[2] and b[4] == height:
                        hist_wat.append(v[11])
                        hist_era.append(b[5])  
                    else:
                        continue
                    
            else:
                continue
        
          
        output= '{}/speed_height_wat_era5_sametimes_scatter_height{}'.format(output_path,height)
        titel='Windspeed for Height: {}m'.format(height) 
        plot_scatter_2D(i_height, hist_wat, hist_era, savefile=output, name_1=name_input_1, name_2=name_input_2, xlabel=name_input_1, ylabel= name_input_2,title=titel,)
        print('plot {} done'.format(i_height))        
        hist_wat=[] 
        hist_era=[]  
        
def plot_wat_lidars_2D(input_wat, heights, name_input_1='lidar1', name_input_2='lidar2',output_path='./plots/testplots'):
    
    hist_lid1=[]
    hist_lid2=[]
    dum1=[]
    dum2=[]
    
    heights=[[555,552],[776,757],[1144,1167],[1365,1371]]
    
    
    for i_height,height in enumerate(heights):
        
        
        for c,v in enumerate(input_wat):
            if c==0: 
                time=v[1]
                dum1.append(v)
            if v[1]==time:
                dum1.append(v)
           
            if v[1]!=time:
                # ~ print(dum1)
                for i_dum,dum in enumerate(dum1):
                    dum2.append(dum[10])
                    
                if ( height[0] in dum2)  and (height[1] in dum2):
                    for a,b in enumerate(dum1):
                        if b[10] == height[0] :
                            hist_lid1.append(b[11])
                        if b[10] == height[1] :
                            hist_lid2.append(b[11])  
                    
                dum1=[]
                dum2=[]
                time=v[1]
            
            if c==len(input_wat)-1:
                for i_dum,dum in enumerate(dum1):
                    dum2.append(dum[10])
                    
                if ( height[0] in dum2)  and (height[1] in dum2):
                    for a,b in enumerate(dum1):
                        if b[10] == height[0] :
                            hist_lid1.append(b[11])
                        if b[10] == height[1] :
                            hist_lid2.append(b[11])  
                    
                dum1=[]
                dum2=[]
                time=v[1]
            
         
        output= '{}/speed_height_wat_era5_sametimes_height{}'.format(output_path,height)
        titel='Lidar1 vs Lidar2' 
        plot_scatter_2D(i_height, hist_lid1, hist_lid2, savefile=output, name_1=name_input_1, name_2=name_input_2, xlabel='lidar1_{}m'.format(height[0]), ylabel= 'lidar2_{}m'.format(height[1]),title=titel,)
        print('plot {} done'.format(i_height))        
        hist_lid1=[] 
        hist_lid2=[]  

def plot_wat_lidars_delta(input_wat, heights,origin_hist1=None, name_input_1='lidar1', name_input_2='lidar2',output_path='./plots/testplots',Nbins=100,Binrange=(-10,10)):
    
    hist_lid1=[]
    hist_lid2=[]
    delta_lid=[]
    dum1=[]
    dum2=[]
    
    heights=[[555,552],[776,757],[1144,1167],[1365,1371]]
    
    
    for i_height,height in enumerate(heights):
        
        
        for c,v in enumerate(input_wat):
            if c==0: 
                time=v[1]
                dum1.append(v)
            if v[1]==time:
                dum1.append(v)
           
            if v[1]!=time:
                # ~ print(dum1)
                for i_dum,dum in enumerate(dum1):
                    dum2.append(dum[10])
                    
                if ( height[0] in dum2)  and (height[1] in dum2):
                    for a,b in enumerate(dum1):
                        if b[10] == height[0] :
                            hist_lid1.append(b[11])
                        if b[10] == height[1] :
                            hist_lid2.append(b[11])  
                    
                dum1=[]
                dum2=[]
                time=v[1]
            
            if c==len(input_wat)-1:
                dum1.append(v)
                for i_dum,dum in enumerate(dum1):
                    dum2.append(dum[10])
                    
                if ( height[0] in dum2)  and (height[1] in dum2):
                    for a,b in enumerate(dum1):
                        if b[10] == height[0] :
                            hist_lid1.append(b[11])
                        if b[10] == height[1] :
                            hist_lid2.append(b[11])  
                    
                dum1=[]
                dum2=[]
                time=v[1]
            
        for a,b in enumerate(hist_lid1):
            delta_lid.append(b-hist_lid2[a])
         
            
        output= '{}/lidar1minuslidar2{}'.format(output_path,height)
        titel='lidar1_{}m_lidar2_{}m'.format(height[0],height[1]) 
        plot_hist(i_height, delta_lid, origin_hist_1=origin_hist1,savefile=output, name_1=name_input_1, xlabel='lidar1-lidar2',nbins=Nbins,binrange=Binrange,title=titel)
        print('plot {} done'.format(i_height))        
        
        hist_lid1=[] 
        hist_lid2=[]  
        delta_lid=[]
                               
def nc_input_to_list(length,date,time,time_seconds,height,speed,direction,vel_u,vel_v,vel_w,met_qc,height_mean,vel_2d):
    
    all_lines=[]
    dum=[]
    for a in range(length):
        dum.append(date[a])
        dum.append(time[a])
        dum.append(time_seconds[a])
        dum.append(height[a])
        dum.append(speed[a])
        dum.append(direction[a])
        dum.append(vel_u[a])
        dum.append(vel_v[a])
        dum.append(vel_w[a])
        dum.append(met_qc[a])
        dum.append(height_mean[a])
        dum.append(vel_2d[a])
        all_lines.append(dum)
        dum=[]
        
    return all_lines

# ~ def calculate_rms(all_lines):
    

    
    
def main():
    
    # ~ #read in wat_data sets
    all_lines_079 = read_wat_data('w_data_079.nc')
    all_lines_0_cut_total = read_wat_data('w_data_0_cut_total.nc')
    
    # ~ #sort w_data sets by height (This has to be done for the delta plots i think...)
    all_lines_079_sort_height = sort_by_index(all_lines_079)
    all_lines_0_cut_total_sort_height = sort_by_index(all_lines_0_cut_total)
    
    
    #read in era5_data
    era5_data = read_wat_era5_data('era5_data_sorted.nc')
    era5_data_sort_height = sort_by_index(era5_data, index=4)
    
    print('Era5 data is read in')
    
    
    # ~ plot_delta_vel_height(all_lines_0, heights_meaned, name_input_1='Wat data', output_path='./plots/cut total singleplots',Nbins=200,Binrange=(-10,10)) 
    # ~ plot_delta_vel_height(all_lines_0_cut_total, heights_meaned, name_input_1='Wat data', output_path='./plots/cut total singleplots',Nbins=200,Binrange=(-10,10)) 

            
    # look at sorted output
    # ~ read_out_data (all_lines_079_sort_height,'a.txt')
    # ~ read_out_data (all_lines_07_sort_height,'b.txt')
    # ~ read_out_data (all_lines_0_sort_height,'c.txt')
    # ~ read_out_data (all_lines_0_cut_jumps_sort_height,'d.txt')
    # ~ read_out_data (all_lines_0_cut_total_sort_height,'e.txt')
    
    #look at output
    # ~ read_out_data (all_lines_079,'a.txt')
    # ~ read_out_data (all_lines_07,'b.txt')
    # ~ read_out_data (all_lines_0,'c.txt')
    # ~ read_out_data (all_lines_0_cut_jumps,'d.txt')
    
    
    
    
    # Plotting era5 vs wattisham
    
    
    
    # ~ plot_watera5_deltarel_vel_height(all_lines_0_cut_total_sort_height, era5_data_sort_height, heights_meaned, name_input_1='(wat-ERA5)/wat',output_path='./plots/cuttotal vs era5 relative',Nbins=50,Binrange=(-1.5,1.5))
    
    # ~ plot_watera5_delta_vel_height(all_lines_0_cut_total_sort_height, era5_data_sort_height, heights_meaned, name_input_1='wat - ERA5',output_path='./plots/cuttotal vs era5',Nbins=100,Binrange=(-20,20))

    
    #compare distributions of wat and era5 wind speeds     
    # ~ plot_watera5_vel_height_2_same_times(all_lines_0_cut_total_sort_height, era5_data_sort_height, heights_meaned, name_input_1='Wat', name_input_2 = 'ERA5', output_path='./plots/era5 and cuts_total_sametimes',Nbins=100,Binrange=(-2,35))
    
    
    
    
    #2D Plots
    # ~ plot_watera5_vel_height_2_same_times_2D(all_lines_0_cut_total_sort_height, era5_data_sort_height, heights_meaned,name_input_1='v_wat', name_input_2='v_era5',output_path='./plots/testplots')
    
    
    
    #Lidar plots
    
    # ~ plot_wat_lidars_2D(all_lines_0_cut_total, heights_meaned, name_input_1='lidar1', name_input_2='lidar2',output_path='./plots/lidar1vs2')
    
    # ~ plot_wat_lidars_delta(all_lines_0_cut_total, heights_meaned,output_path='./plots/lidar1vs2',Nbins=100,Binrange=(-10,10))
    
    
    
    
    # ~ #three speed areas, wat vs era5
    # ~ l1=[]
    # ~ l2=[]
    # ~ l3=[]
    
    # ~ for c,v in enumerate(all_lines_0_cut_total_sort_height):
        # ~ if v[11] < 8:
            # ~ l1.append(v)
        # ~ if v[11] > 8 and v[11] < 12:
            # ~ l2.append(v)
        # ~ if v[11] > 12:
            # ~ l3.append(v)
    
    # ~ plot_watera5_delta_vel_height(l1, era5_data_sort_height, heights_meaned, name_input_1='wat-era5',output_path='./plots/testplots/l1')
    # ~ plot_watera5_delta_vel_height(l2, era5_data_sort_height, heights_meaned, name_input_1='wat-era5',output_path='./plots/testplots/l2')
    # ~ plot_watera5_delta_vel_height(l3, era5_data_sort_height, heights_meaned, name_input_1='wat-era5',output_path='./plots/testplots/l3')
    
    

    
    #Hists to justify the cuts vel=4 and direction=10
    # ~ plot_delta_vel_height(all_lines_0, heights_meaned, name_input_1='Wat data', output_path='./plots/cut total singleplots',Nbins=200,Binrange=(-10,10)) 
    # ~ plot_delta_dir_height(all_lines_0, heights_meaned, name_input_1='Wat data', output_path='./plots/cut total singleplots',Nbins=60,Binrange=(-30,30))
    # ~ plot_delta_dir_time(all_lines_0_cut_total_sort_height, times, name_input_1='cuts_total')
    # ~ plot_delta_vel_time(all_lines_0_cut_total_sort_height, times, name_input_1='cuts_total') 
     
    
    # qc_079 vs qc_0
    # ~ plot_delta_vel_time_2(all_lines_079_sort_height,all_lines_0_sort_height, times, name_input_1='079', name_input_2='0', output_path='./plots/qc_079 vs qc_0')
    # ~ plot_delta_vel_height_2(all_lines_079,all_lines_0, heights_meaned, name_input_1='079', name_input_2='0', output_path='./plots/qc_079 vs qc_0')
    # ~ plot_delta_dir_time_2(all_lines_079_sort_height,all_lines_0_sort_height, times, name_input_1='079', name_input_2='0', output_path='./plots/qc_079 vs qc_0')
    # ~ plot_delta_dir_height_2(all_lines_079,all_lines_0, heights_meaned, name_input_1='079', name_input_2='0', output_path='./plots/qc_079 vs qc_0')
    
    # ~ #qc_0_cuts_jumps vs qc_0
    # ~ plot_delta_vel_time_2(all_lines_0_cut_jumps_sort_height,all_lines_0_sort_height, times, name_input_1='qc_0_cut_jumps', name_input_2='qc_0', output_path='./plots/qc_0_cuts_jumps vs qc_0')
    # ~ plot_delta_vel_height_2(all_lines_0_cut_jumps,all_lines_0, heights_meaned, name_input_1='qc_0_cut_jumps', name_input_2='qc_0', output_path='./plots/qc_0_cuts_jumps vs qc_0')
    # ~ plot_delta_dir_time_2(all_lines_0_cut_jumps_sort_height,all_lines_0_sort_height, times, name_input_1='qc_0_cut_jumps', name_input_2='qc_0', output_path='./plots/qc_0_cuts_jumps vs qc_0')
    # ~ plot_delta_dir_height_2(all_lines_0_cut_jumps,all_lines_0, heights_meaned, name_input_1='qc_0_cut_jumps', name_input_2='qc_0', output_path='./plots/qc_0_cuts_jumps vs qc_0')
    
    # ~ #qc_07 vs qc_0
    # ~ plot_delta_vel_time_2(all_lines_07_sort_height,all_lines_0_sort_height, times, name_input_1='07', name_input_2='0', output_path='./plots/qc_07 vs qc_0')
    # ~ plot_delta_vel_height_2(all_lines_07,all_lines_0, heights_meaned, name_input_1='07', name_input_2='0', output_path='./plots/qc_07 vs qc_0')
    # ~ plot_delta_dir_time_2(all_lines_07_sort_height,all_lines_0_sort_height, times, name_input_1='07', name_input_2='0', output_path='./plots/qc_07 vs qc_0')
    # ~ plot_delta_dir_height_2(all_lines_07,all_lines_0, heights_meaned, name_input_1='07', name_input_2='0', output_path='./plots/qc_07 vs qc_0')
    
    # ~ #qc_cuts_total vs qc_0
    # ~ plot_delta_vel_time_2(all_lines_0_cut_total_sort_height,all_lines_0_sort_height, times, name_input_1='cuts_total', name_input_2='uncut', output_path='./plots/cuts_total vs qc_0',Nbins=200,Binrange=(-10,10))
    # ~ plot_delta_dir_time_2(all_lines_0_cut_total_sort_height,all_lines_0_sort_height, times, name_input_1='cuts_total', name_input_2='uncut', output_path='./plots/cuts_total vs qc_0',Nbins=200,Binrange=(-10,10))
    plot_delta_vel_height_2(all_lines_0_cut_total,all_lines_079, heights_meaned, name_input_1='wat with cuts', name_input_2='wat uncut', output_path='./plots/cuts_total vs qc_079',Nbins=200,Binrange=(-10,10))
    plot_delta_dir_height_2(all_lines_0_cut_total,all_lines_079, heights_meaned, name_input_1='wat with cuts', name_input_2='wat uncut', output_path='./plots/cuts_total vs qc_079',Nbins=60,Binrange=(-30,30))
    
     
    
    
    
    
           
if __name__ == "__main__":
    main()
