'''
Created ab4810 on 27/01/2020

Python script to generate raypath mask for AARM derived datasets.

This is parallelized code for computing mask grids from the PP files.
Not written in nested functions as cannot be pickled for passage to/from cores.

Dataset in form:
evlat, evlon, evdep, stlat, stlon, arc_d, pp_dep, pp_lat, pp_lon, pp_fzhw

Here we want to calculate the mask grid
specify file-in_loc, file-in filename, output mask_grid_file, lat/lon bounds and spacial sampling interval
fzmf = 'Fresnel zone multiplication factor' Used to broaden the footprint of each fresnel zone based mask.
fzgs = 'Fresnel zone grid size' A static value added to region "captured" by ray. Tomography is gridded in constant velocity cells of minimum size ~fzgs
        This could be modified to be depth adaptive based on average tomographic grid cell density in region etc.



cores=4
home = '/Users/ab4810'
file_in_loc='/Users/ab4810/Google_Drive/GITHUB_AB/3D_model_info'
PP_FZ_file_in='AFRP20_phase_SUMMARY_PP_FZ.txt'
mask_grid_file='AFRP20_MASK'
lon_min=-24.0
lon_max=64.0
lat_min=-44.0
lat_max=44.0
grid_d_lat_lon=0.25
fzmf=1.0
fzgs=40.0

file_in_loc='/Users/ab4810/Google_Drive/GITHUB_AB/3D_model_info'
PP_FZ_file_in='ATS_50_RAYS_PIERCE.txt'
mask_grid_file='TEST_ATS_50_MASK'
lon_min=-24.0
lon_max=64.0
lat_min=-44.0
lat_max=44.0
grid_d_lat_lon=0.25
fzmf=1.0
fzgs=40.0


'''
### Importing various python libraries
# numpy is a useful toolkit for scientific computations
import numpy as np
import sys,glob
import os.path
import math
import time
import concurrent.futures

print('Imports complete.....')

cores=4
home = '/raid1/ab2568'
file_in_loc=home+'/GITHUB_AB/3D_model_info'
PP_FZ_file_in='AFRP20_phase_SUMMARY_PP_FZ.txt'
mask_grid_file='AFRP20_MASK_1.0'
lon_min=-24.0
lon_max=64.0
lat_min=-44.0
lat_max=44.0
grid_d_lat_lon=1
fzmf=1.0
fzgs=40.0

def haversine(lat1, long1, lats2, longs2, depth):
    """
    Input depth in km
    Calculate the distance between two points in earth output in km
    """
    d=[]

    for i in range(len(lats2)):
        lat2=lats2[i]
        long2=longs2[i]
        radius = 6371.e3 - depth*1.e3  # m
        dLat = math.radians(lat2 - lat1)
        dLong = math.radians(long2 - long1)

        a = (math.sin(dLat / 2.0) ** 2.0 +
             math.cos(math.radians(lat1)) * math.cos(math.radians(lat2)) * math.sin(dLong / 2.0) ** 2.0)
        c = 2.0 * math.atan2(math.sqrt(a), math.sqrt(1.0 - a))
        d.append(radius * c/1000.0) # div 1000 for kms
    return np.array(d)



######## Check input parameters ###############################
if not isinstance(file_in_loc, str):
    print('Bad specification of file_in_loc - must be string')
    sys.exit()
if not isinstance(PP_FZ_file_in, str):
    print('Bad specification of filename - must be string')
    sys.exit()
if not os.path.exists(file_in_loc+'/'+PP_FZ_file_in):
    print('Bad specification of file_in_loc + PP_FZ_file_in - does not exist!')
    sys.exit()
if not isinstance(mask_grid_file, str):
    print('Bad specification of out_file - must be string')
    sys.exit()
    
##############################################################

# First make the outfile
dep_mask_filename=file_in_loc + '/' + mask_grid_file+'.txt'
if os.path.exists(dep_mask_filename):
    print('Bad specification of out_file :' + str(dep_mask_filename))
    try:
        print('deleting old out-file: '+str(dep_mask_filename)+'...')
        os.remove(dep_mask_filename)
    except:
        print('Failed to delete '+str(dep_mask_filename))
        sys.exit()
dep_mask_file = open(dep_mask_filename, 'w')
dep_mask_file.write('{0:2s} {1:3s} {2:4s} {3:5f} {4:6f} {5:7f} {6:8f} {7:9f} {8:10f} {9:11f} \n'.format(file_in_loc, PP_FZ_file_in, mask_grid_file,lon_min,lon_max,lat_min,lat_max,grid_d_lat_lon,fzmf,fzgs))
###########################################################

print('Starting grid calcs...')

data_file = open(file_in_loc+'/'+PP_FZ_file_in,'r')
d_vector=data_file.read().strip().split()
d_vector = [float(x) for x in d_vector]
data_list=np.reshape(d_vector,(-1,10)) # 10 is number of columns in infile

# Make an array with interogation depths, lats and lons
dep_ints=np.arange(0,2891,50)
g_lats=np.arange(lat_min,lat_max+grid_d_lat_lon,grid_d_lat_lon)
g_lons=np.arange(lon_min,lon_max+grid_d_lat_lon,grid_d_lat_lon)

grid_points=[]
lat_points=[]
lon_points=[]
for i in range(len(g_lats)):
    for j in range(len(g_lons)):
        lat_points.append(g_lats[i])
        lon_points.append(g_lons[j])
lon_points=np.array(lon_points)
lat_points=np.array(lat_points)


#Loop over the dep interval array
for dep in dep_ints:
    dep_start = time.perf_counter()
    # Make dep array to paste with each mask
    dep_points=dep*np.ones(len(lat_points))
    
    # Find rows of data_list that correspond to correct depth.
    dep_index=np.where(data_list[:,6]==dep)[0]
    num_pp=len(dep_index)
    print('Found '+str(num_pp)+' pierce points at '+str(dep)+'km depth..')
    # For each selected row of data_list
    if num_pp>1:
        # for ind in dep_index:
        def get_grids(dep_index):
            global data_list, lat_points, lon_points, dep, fzmf, fzgs
            # use pierce point and Fresnel zone
            pp_lat=data_list[dep_index,7]
            pp_lon=data_list[dep_index,8]
            fzhw  =data_list[dep_index,9]
            # Calculate distances to all points on the grid
            grid_to_pp_dist=haversine(pp_lat, pp_lon, lat_points, lon_points, dep)
            # Make a mask grid "on" everywhere - i.e. equal to one
            mask_grid=np.ones(len(lat_points))
            # Convert all points in mask grid less than specified fresnel zone buffer distance to zero - mask off.
            mask_grid[grid_to_pp_dist<=(fzhw*fzmf+fzgs)]=0
            # Add these mask grids to the array dep_mask_grid for each pp
            return mask_grid

        # This method adds as we go along so should run into less memory trouble.
        dep_mask_grid_out=np.zeros(len(lat_points))
        with concurrent.futures.ProcessPoolExecutor(max_workers=cores) as executor:
            results = executor.map(get_grids, dep_index)
            # Loop over each result and add to developing depth mask
            for mask_grid in results:
                add_grid=np.array(mask_grid)
                # want to sum each individual mask value along same dimension
                dep_mask_grid_out=dep_mask_grid_out+add_grid
        
        dep_mask_grid_out[dep_mask_grid_out<num_pp]=0
        dep_mask_grid_out[dep_mask_grid_out==num_pp]=1
        
        
        print('Writing to file: depth: '+str(dep)+'km')
        for i in range(len(lat_points)):
            dep_mask_file.write('{0:2f} {1:3f} {2:4f} {3:5f}\n'.format(dep_points[i], lat_points[i], lon_points[i], dep_mask_grid_out[i]))
        dep_finish= time.perf_counter()

        print(f'Finished in {round(dep_finish-dep_start, 2)} seconds(s)\n')
dep_mask_file.close()
