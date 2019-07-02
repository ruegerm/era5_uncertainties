import os
from netCDF4 import Dataset, MFDataset
import numpy as np

from config_wat_analyse import wat_data_dir, output_text_file, output_nc_file


req_altitudes = [184, 256, 328, 401,473, 545, 617, 690, 762, 834, 907, 979, 1051, 1124, 1196, 1268, 1341, 1413, 1485, 1557, 1630, 1702]
                       
req_altitudes2 = [184, 190, 256, 265, 328, 340, 347, 401, 415, 473, 490, 545, 552, 565, 617, 640, 690, 715, 757, 762, 790, 
834, 865, 907, 940, 962, 979, 1014, 1051, 1089, 1124, 1164, 1167, 1196, 1239, 1268, 1314, 1341, 1371, 1389,
 1413, 1464, 1485, 1539, 1557, 1576, 1614, 1630, 1689, 1702,1764]

#Lidar 1 operated at two different heights according to the data, therefore the mean of the corresponding height pairs is taken. The data is evaluated at these heights using linear interpolation
heights_lidar1a=['184', '256', '328', '401', '473', '545', '617', '690', '762', '834', '907', '979', '1051', '1124', '1196', '1268', '1341', '1413', '1485', '1557', '1630']
heights_lidar1b=['190', '265', '340', '415', '490', '565', '640', '715', '790', '865', '940', '1014', '1089', '1164', '1239', '1314', '1389', '1464', '1539', '1614', '1689']

heights_lidar1_mean = [187, 260, 334, 408, 481, 555, 628, 702, 776, 849, 923, 996, 1070, 1144, 1217, 1291, 1365, 1438, 1512, 1585, 1659]
heights_lidar2 = [347, 552, 757, 962, 1167, 1371, 1576]
heights_all_meaned= [187, 260, 334, 347, 408, 481, 552, 555, 628, 702, 757, 776, 849, 923, 962, 996, 1070, 1144, 1167, 1217, 1291, 1365, 1371, 1438, 1512, 1576, 1585, 1659]

# Wat data indices: 0 = date, 1 = time, 2= time_secs, 3 = height, 4 = total windspeed, 5 = wind direction, 6 = vel_x, 7 = vel_y, 8 = vel_z, 9 = metqc, 10 = mean_heights, 11 = vel_2d
# Era5 data indices: 0 = hours_since_1900, 1 = date, 2 = time, 3 = altitude_sea, 4 = altitude_ground = wat_height_mean, 5 = wind_speed

def convert_int_to_time(time_int):
    # ~ # convert '080000' time value to time_in_seconds_since_midnight
    a = [time_int[0:2],time_int[2:4],time_int[4:6]]
    t_seconds = int(a[0])*3600 + int(a[1])*60 +int(a[2])
    return t_seconds
   
def read_watthisam_file(input_filename, output_filename):
    # this function reads in an original watthisam_data_file and writes the relevant information to an output file.txt
    
    with open (input_filename, 'r') as recent_file:  
        all_lines=[]
        
        # creating a list of every measurement in input_filename
        for line in recent_file:    
            all_lines.append(line.split())

    # analyzing the list and only read out the interesting stuff (date, time, wind data, met_qc etc.)
    with open (output_filename, 'a') as output:
        for c,v in enumerate(all_lines):
           
            
            if c==0:
                continue
            if 'Wattisham' in v:
                date=all_lines[c+3][0]+all_lines[c+3][1]+all_lines[c+3][2]
                time=all_lines[c+3][3]+all_lines[c+3][4]+all_lines[c+3][5]
                
                continue
            if 'HT' in v:
                continue
            if len(v)>11:
                l_dummy=[date,time,convert_int_to_time(time)]+v[0:7]
                
                for i,j in enumerate(l_dummy):
                    if i==3:    #change height to unit meters
                        a=int(float(j)*1000)
                        output.write('%s\t' %a)
                    else:
                        output.write('%s\t' %j)
                    
                output.write('\n')
                                    
                    
                    
def wat_data_to_netcdf(input_file_name,output_file):
    # ~ # this function creates a netcdf file from an input text file
    
    print('Now creating netcdf file...') 
    dates = []
    times = []
    times_seconds =[]
    heights = []
    speeds= []
    directions= [] 
    vel_us = [] 
    vel_vs = [] 
    vel_ws = [] 
    met_qcs = []
    vel_2ds = []    #2d projection of velocity vector
    heights_mean = []
    
    # ~ read in the data
    with open (input_file_name, 'r') as recent_file:  
        all_lines=[]
        
        for c,line in enumerate(recent_file):    
            all_lines.append(line.split())
    
       
    # ~ open netcdf Dataset and create variables            
    nc_out = Dataset(output_file, "w", format="NETCDF3_64BIT_OFFSET")
    nc_out.createDimension("length", len(all_lines))
    
    date = nc_out.createVariable("date", "i4", ("length",))
    time = nc_out.createVariable("time", "i4", ("length",))
    time_seconds = nc_out.createVariable("time_seconds", "i4", ("length",))
    height = nc_out.createVariable("height", "i4", ("length",))
    speed = nc_out.createVariable("speed", "f4", ("length",))
    direction = nc_out.createVariable("direction", "i4", ("length",))
    vel_u = nc_out.createVariable("vel_u", "f4", ("length",))
    vel_v = nc_out.createVariable("vel_v", "f4", ("length",))
    vel_w = nc_out.createVariable("vel_w", "f4", ("length",))
    met_qc = nc_out.createVariable("met_qc", "i4", ("length",))
    vel_2d = nc_out.createVariable("vel_2d", "f4", ("length",))
    height_mean = nc_out.createVariable("height_mean", "i4", ("length",))
    
    # split data in all_lines into different variables 
    for c,v in enumerate(all_lines):
        dates.append(v[0])
        
        times.append(v[1])
        times_seconds.append(v[2])
        heights.append(v[3])
        speeds.append(v[4])
        directions.append(v[5])
        vel_us.append(v[6])
        vel_vs.append(v[7])
        vel_ws.append(v[8])
        met_qcs.append(v[9])
        heights_mean.append(v[10])
        vel_2ds.append(v[11])
        
    # put data into netcdf variables     
    date[:]=dates
    time[:]=times
    time_seconds[:]=times_seconds
    height[:]=heights
    speed[:]= speeds
    direction[:]= directions
    vel_u[:]= vel_us
    vel_v[:]= vel_vs
    vel_w[:]= vel_ws
    met_qc[:]= met_qcs
    vel_2d[:]= vel_2ds
    height_mean[:]= heights_mean
    
    nc_out.close()
    print('Netcdf file is created') 

 
    
def correct_full_hours(all_lines):
    #correcting the measurements with time values at full hours by setting their time value to somewhere in the hour before so that the analysis code can identify them correctly
    all_lines_new=[]
    
    for c,v in enumerate(all_lines):
        
        if c==0:
            all_lines_new.append(v)
            continue
        
        if all_lines[c][0]!=all_lines[c-1][0] and int(all_lines[c][2])==0: #if date not the same and time == 0
            dum1=all_lines[c-1][0]
            all_lines_new.append(v)
            all_lines_new[c][0]=dum1
            all_lines_new[c][1]=235959
            all_lines_new[c][2]=convert_int_to_time('235959')
            continue
            
        if all_lines[c][0]==all_lines[c-1][0] and int(all_lines[c][2])/1000==0: #if date the same and time a full hour
            dum1=all_lines[c-1][0]
            all_lines_new.append(v)
            all_lines_new[c][0]=dum1
            all_lines_new[c][1]='{}5959'.format((int(all_lines_new[c][2][0:2])-1))
            all_lines_new[c][2]=convert_int_to_time('235959')
            continue
            
        all_lines_new.append(v)
    
    return all_lines_new
    
    
    
def extend_data(all_lines):
    #adding the mean height value of lidar1 measurements to all_lines and the 2D windspeed projection
    
    print('Now extending data...')
    all_lines_new=all_lines
    for c,v in enumerate(all_lines):
        
        if v[3] in (heights_lidar1a):
            all_lines_new[c].append(heights_lidar1_mean[heights_lidar1a.index(v[3])])
            all_lines_new[c].append( round((float(v[6])**2+float(v[7])**2)**.5,2))
            continue
        
        if v[3] in (heights_lidar1b):
            all_lines_new[c].append(heights_lidar1_mean[heights_lidar1b.index(v[3])])
            all_lines_new[c].append( round((float(v[6])**2+float(v[7])**2)**.5,2))
            continue
        
        all_lines_new[c].append(v[3])
        all_lines_new[c].append( round((float(v[6])**2+float(v[7])**2)**.5,2))
        
        
    print('Data is extended')
    return all_lines_new
    
def read_in_data(data_file):
    # reading in lines from a txt file into a list
    
    
    with open (data_file, 'r') as recent_file:  
        all_lines=[]
        for line in recent_file:
            all_lines.append(line.split())
    print('Data is read in')
      
    return all_lines
        

def cut_on_height(all_lines,height):
    #cut a wat_data_list on values above a given height
    
    print('Now cutting the height with a knife...')
    all_lines_new=[]
    for c,v in enumerate(all_lines):
        
        if int(v[3])<height: 
            
            all_lines_new.append(v)
    print('Cut on height has been performed')
    return all_lines_new

def read_out_data(input_list,output_file):
    # ~ #read input_list into output_file 
    
    with open(output_file, 'w') as f:
        
        for a in input_list:
            for b in a:
                f.write('%s \t' %b)
            f.write('\n')
    print('Data is read out.')

def sort_index_insideatime(all_lines, index):
    #sorting a wat_data_list by a given index of the list inside one measurement time
    
    print('Sorting the data...')
    l_dummy=[]
    all_lines_new=[]
    
    for c,v in enumerate(all_lines):
        
        if c==0:
            l_dummy.append(v)
            continue
        
        if v[0]==all_lines[c-1][0]: #if date the same 
            
            if ( abs(int(v[2])-int(all_lines[c-1][2])) ) < 600 : #if time the same
                l_dummy.append(v)
            if ( abs(int(v[2])-int(all_lines[c-1][2])) ) > 600 : #if time not the same
                
                l_dummy.sort(key=lambda x:int(x[index-1]))
                all_lines_new.extend(l_dummy)
                
                l_dummy=[]
                l_dummy.append(v)
                continue
        if v[0]!=all_lines[c-1][0]: #if date not the same 
            
            l_dummy.sort(key=lambda x:int(x[index-1]))
            all_lines_new.extend(l_dummy)
            
            l_dummy=[]
            l_dummy.append(v)
            continue
    print('How i love sorted data! This one is totally sorted')          
    return all_lines_new
        
def sort_by_index(all_lines, index=10):
    #sorting the data by any column index, standard is the height column
    
    print('Sorting the data by height...')
    all_lines_new=[]
    for c,v in enumerate(all_lines):
        all_lines_new.append(v)
    all_lines_new.sort(key=lambda x:int(x[index]))
    
    return all_lines_new

def sort_ind1_ind2_ind3(all_lines, index1=10, index2=1, index3=0):
    #sorting the data by index1, then index2, then index3, standard is index1=height, index2=hour, index3=date
    print('Sorting the data by height...')
    all_lines_new=[]
    for c,v in enumerate(all_lines):
        all_lines_new.append(v)
    all_lines_new.sort(key=lambda x:int(x[index1]))
    all_lines_new.sort(key=lambda x:int(x[index2]))
    all_lines_new.sort(key=lambda x:int(x[index3]))
    
    return all_lines_new    
        
def mean_wat_vel(all_lines):
   # for all data get the mean values for all measurements inside one hour
    print('Now merging times...')
    all_lines_new=[]
    
    # typecasting the all_lines list
    for c,v in enumerate(all_lines):
        all_lines[c][4]=float(v[4])
        all_lines[c][6]=float(v[6])
        all_lines[c][7]=float(v[7])
        all_lines[c][8]=float(v[8])
        all_lines[c][11]=float(v[11])
        
    x=30
    all_lines_new=[]
    l_dummy=[]
    for c,v in enumerate(all_lines):
       
        if x==30:
            l_dummy.append(v)
            x=int(v[2])/3600
            continue
            
        if v[0]==all_lines[c-1][0]: #if date the same 
            
            if ( x== int(v[2])/3600 ): #if time the same
                
                l_dummy.append(v)
                
            if ( x!= int(v[2])/3600 ): #if time not the same
                
                l_dummy=merge_times(l_dummy,x)
                all_lines_new.extend(l_dummy)
                l_dummy=[]
                l_dummy.append(v)
                x=int(v[2])/3600
                continue
                
        if v[0]!=all_lines[c-1][0]: #if date not the same 
            
            
            l_dummy=merge_times(l_dummy,x)
            all_lines_new.extend(l_dummy)
            l_dummy=[]
            l_dummy.append(v)
            x=30
            continue 
    
    # "typecasting" the data
    for c,v in enumerate(all_lines_new):
        v[0]= int(v[0])
        v[1]= int(v[1])
        v[2]= int(v[2])
        v[3]= int(v[3])
        v[4]= float(v[4])
        v[5]= int(v[5])
        v[6]= float(v[6])
        v[7]= float(v[7])
        v[8]= float(v[8])
        v[9]= int(v[9])
        v[10]= int(v[10])
        v[11]= float(v[11])
    print('Times merged')
    return all_lines_new

def merge_times(array,x):
    # calculate the mean values of a list which contains all data for one hour
    
    l_dum=[]
    l_dum2=[i for i in range(12)] 
    l_return= []
    l_heights=[]
    
    # fill l_heights with all heights in array
    for c,v in enumerate(array):
        if v[10] not in l_heights:
            l_heights.append(v[10])
        
    for c,v in enumerate(l_heights):
        
        for a,b in enumerate(array):
            if int(b[10])==int(v):
                
                l_dum.append(b)
        #calculate the new mean values for the hour  
        l_dum2[0]=l_dum[0][0]
        l_dum2[1]=x
        l_dum2[2]=x*3600
        l_dum2[3]=l_dum[0][3]
        l_dum2[4]=add_matrix_col(l_dum,4)
        l_dum2[5]=l_dum[0][5]
        l_dum2[6]=add_matrix_col(l_dum,6)
        l_dum2[7]=add_matrix_col(l_dum,7)
        l_dum2[8]=add_matrix_col(l_dum,8)
        l_dum2[9]=l_dum[0][9]
        l_dum2[10]=l_dum[0][10]
        l_dum2[11]=add_matrix_col(l_dum,11)
        
        l_return.append(l_dum2)
        l_dum=[]
        l_dum2=[i for i in range(12)]    
        
    return l_return
        
def add_matrix_col(matrix,i_col):
    # Adding a column in a given matrix
    
    sum_col=0
    n_values=0
    for c,v in enumerate(matrix):
        # check for bad values
        if v[i_col]<10000:
            sum_col += v[i_col]
            n_values+=1
    if n_values==0:
        return 999999
    return round(sum_col/n_values,2)

def cut_by_metqc(all_lines,qc1=None,qc2=None):
    # cut all values with quality met_qc equal to qc1 and/or qc2
    
    print('Starting to clean the data...')
    all_lines_new=[]
    
    for c,v in enumerate(all_lines):
        
        if qc1 != None:
            
            if int(v[9]) == qc1:
                continue
       
        if qc2 != None:
           
           if int(v[9]) == qc2:
               continue
       
        all_lines_new.append(v)
       
    print('Data is so clean, you can eat from it')               
    return all_lines_new

def jumps_closebadmetqc(all_lines,jump_vel,jump_dir):
    #cut all values with Delta = jump_vel to a neighbour value and close to bad met_qc values
    
    print('Now taking jumps away...')
    all_lines_new=[]
    all_lines_new2=[]
    for c,v in enumerate(all_lines):
        
        if c==0:
            all_lines_new.append(v)
            continue
        
        if c==len(all_lines)-1:
            if (abs(v[4]-all_lines[c-1][4]) > jump_vel or abs(get_distance_directions(v[5],all_lines[c-1][5])) >jump_dir )  and (all_lines[c-1][9] !=0):
                break
            all_lines_new.append(v)
            break
        
        if v[0]==all_lines[c-1][0]: #if date the same 
            
            if  v[2]== all_lines[c-1]  : #if time the same
                
                if (abs(v[4]-all_lines[c-1][4]) > jump_vel and all_lines[c-1][4]!=999999) and all_lines[c+1][9] !=0:
                    continue
                if (abs(v[4]-all_lines[c+1][4]) > jump_vel and all_lines[c+1][4]!=999999) and all_lines[c-1][9] !=0:
                    continue
                if (abs(get_distance_directions(v[5],all_lines[c-1][5])) > jump_dir and all_lines[c-1][4]!=999999  ) > jump_dir and all_lines[c+1][9] !=0:
                    continue 
                if (abs(get_distance_directions(v[5],all_lines[c+1][5])) > jump_dir and all_lines[c+1][4]!=999999  ) > jump_dir and all_lines[c-1][9] !=0:
                    continue
        
                all_lines_new.append(v)
                
            if v[2]!= all_lines[c-1] : #if time not the same
                all_lines_new.append(v)
                continue
                
        if v[0]!=all_lines[c-1][0]: #if date not the same 
            all_lines_new.append(v)
            continue
    
    return all_lines_new


def jumps_bothneighbor_heights(all_lines, jump_vel, jump_dir):
     # cut values with jumps to both neighbours in windspeed and direction
    all_lines_new=[]
    for c,v in enumerate(all_lines):
        if c==0:
            all_lines_new.append(v)
            continue
        
        if c==len(all_lines)-1:
            all_lines_new.append(v)
            break
       
        if (abs(v[4]-all_lines[c-1][4]) > jump_vel and abs(v[4]-all_lines[c+1][4]) > jump_vel) and v[9] ==0:
            
            continue
        
        if (abs(get_distance_directions(v[5],all_lines[c-1][5])) > jump_dir and abs(get_distance_directions(v[5],all_lines[c+1][5])) > jump_dir  ) and v[9]==0:
            
            continue 
        
        all_lines_new.append(v)
    
    return all_lines_new
    
    
    
def cut_999999_values(all_lines):
    # cut measurements which contain any corrupted value
    all_lines_new=[]
    for c,v in enumerate(all_lines):
        
        if 999999 in v:
            continue
        all_lines_new.append(v) 
    return all_lines_new
    
def get_distance_directions(a,b):
    # calculate the absolute deviation of two directions given in degrees 
    
    if (abs(int(a)-int(b))>180):
        c= 360-(int(a)-int(b))
    else:
        c=int(a)-int(b)
    return c
         
def omit_doubles(all_lines):
    # get rid of any double data for one height 
    
    all_lines_new=[]
    for c,v in enumerate(all_lines):         
        if c== 0 or c== len(all_lines)-1:
            all_lines_new.append(v)
            continue
        if v[10]== all_lines[c-1][10]:
            continue
        else:
            all_lines_new.append(v)
    return all_lines_new

def arrange_data(all_lines):
   # arrange the data in a convenient way
   
    data_cut_height = cut_on_height(all_lines,1750)
    data_extended = extend_data(data_cut_height)
    data_sorted_by_time = sort_index_insideatime(data_extended,11)
    data_corrected = correct_full_hours(data_sorted_by_time)
    data_meaned = mean_wat_vel(data_corrected)

    return data_meaned

def filter_data(all_lines):
    # filter data for suspicious data points
    
    data_0 = cut_by_metqc(all_lines,7,9)
    data_cut_jumps1 = jumps_closebadmetqc(data_0,4,10)
    data_cut_jumps2 = jumps_bothneighbor_heights(data_cut_jumps1,4,10)
    data_0_cut_all = cut_999999_values(data_cut_jumps2)
    
    return data_0_cut_all

def create_raw_wat_file():
    # Read in all wattisham files and write all the data into one big text_file
     
    #collect the data files
    l_files = []
    for root, dirs, filenames in os.walk(wat_data_dir):   # read in datafile names
        for f in filenames:
            l_files.append(os.path.join(root, f))
    
    #sorting the files by name, i.e. date
    l_files.sort()    
    
    #create the raw_output_text_file
    with open(output_text_file, 'w') as f:
        pass
    
    #loop over all files to read the data into a text file
    
    for file in l_files:
        read_watthisam_file(file,output_text_file)
        
    print('Raw Wattisham file is created.')   
    
         
        
def main():
    
    
    
    # ~ #Create the raw wattisham text file
    # ~ create_raw_wat_file()
    
    
    data_raw= read_in_data(output_text_file)
    
    #arrange the data in a convenient way (cut on height, take mean of hours, sort etc.)
    data_arranged = arrange_data(data_raw)
    
    # ~ #read out unfiltered wat_data
    read_out_data(data_arranged, 'w_data_079.txt')
    wat_data_to_netcdf('w_data_079.txt','w_data_079.nc')
    
    
    #filter data for suspicious values (discard bad metqc and 999999 data, cut on velocity(4m/s) and direction(10m/s) (
    data_filtered = filter_data(data_arranged)
    
    #read out filtered wat_data
    read_out_data(data_filtered, 'w_data_0_cut_total.txt')
    wat_data_to_netcdf('w_data_0_cut_total.txt','w_data_0_cut_total.nc')
    
    
         
if __name__ == "__main__":
    main()
        
        
         

        
        
   




        
    
